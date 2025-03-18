/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 *  @file   DDMarlinPandora/src/DDPfoCreator.cc
 *
 *  @brief  Implementation of the pfo creator class.
 *
 *  $Log: $
 */

#include "DDPfoCreator.h"
#include "CalorimeterHitType.h"

#include "Api/PandoraApi.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/Vertex.h"
#include "edm4hep/VertexCollection.h"

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>

DDPfoCreator::DDPfoCreator(const Settings& settings, const pandora::Pandora* const pPandora, IMessageSvc* msgSvc)
    : m_settings(settings), m_pandora(*pPandora), m_msgSvc(msgSvc) {}

//------------------------------------------------------------------------------------------------------------------------------------------

DDPfoCreator::~DDPfoCreator() {}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDPfoCreator::CreateParticleFlowObjects(
    edm4hep::ClusterCollection&                     pClusterCollection,
    edm4hep::ReconstructedParticleCollection&       pReconstructedParticleCollection,
    edm4hep::VertexCollection&                      pStartVertexCollection,
    std::vector<const edm4hep::MCParticleCollection*>     MCParticleCollections,
    std::vector<const edm4hep::TrackCollection*>          trackCollections,
    std::vector<const edm4hep::CalorimeterHitCollection*> eCalCollections,
    std::vector<const edm4hep::CalorimeterHitCollection*> hCalCollections,
    std::vector<const edm4hep::CalorimeterHitCollection*> mCalCollections,
    std::vector<const edm4hep::CalorimeterHitCollection*> lCalCollections,
    std::vector<const edm4hep::CalorimeterHitCollection*> lhCalCollections) {

  const pandora::PfoList* pPandoraPfoList = NULL;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(m_pandora, pPandoraPfoList));

  pandora::StringVector subDetectorNames;
  this->InitialiseSubDetectorNames(subDetectorNames);

  std::unordered_map<int, const edm4hep::CalorimeterHitCollection*> calorimeterCollectionMap;
  auto processCollections = [&](const std::vector<const edm4hep::CalorimeterHitCollection*>& collections) {
    for (auto* collection : collections) { calorimeterCollectionMap[collection->getID()] = collection; }
  };
  processCollections(eCalCollections);
  processCollections(hCalCollections);
  processCollections(mCalCollections);
  processCollections(lCalCollections);
  processCollections(lhCalCollections);
  std::unordered_map<int, const edm4hep::TrackCollection*> trackCollectionMap;
  for (auto* collection : trackCollections) { trackCollectionMap[collection->getID()] = collection; }

  //pClusterCollection.parameters().setValues("ClusterSubdetectorNames", subDetectorNames);

  // Create lcio "reconstructed particles" from the pandora "particle flow objects"
  for (pandora::PfoList::const_iterator pIter = pPandoraPfoList->begin(), pIterEnd = pPandoraPfoList->end();
       pIter != pIterEnd; ++pIter) {
    const pandora::ParticleFlowObject* const pPandoraPfo(*pIter);
    edm4hep::MutableReconstructedParticle  pReconstructedParticle0 = pReconstructedParticleCollection.create();
    edm4hep::MutableReconstructedParticle* pReconstructedParticle  = &pReconstructedParticle0;

    const bool                  hasTrack(!pPandoraPfo->GetTrackList().empty());
    const pandora::ClusterList& clusterList(pPandoraPfo->GetClusterList());

    float                    clustersTotalEnergy(0.f);
    pandora::CartesianVector referencePoint(0.f, 0.f, 0.f), clustersWeightedPosition(0.f, 0.f, 0.f);
    for (pandora::ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end();
         cIter != cIterEnd; ++cIter) {
      const pandora::Cluster* const pPandoraCluster(*cIter);
      pandora::CaloHitList          pandoraCaloHitList;
      pPandoraCluster->GetOrderedCaloHitList().FillCaloHitList(pandoraCaloHitList);
      pandoraCaloHitList.insert(pandoraCaloHitList.end(), pPandoraCluster->GetIsolatedCaloHitList().begin(),
                                pPandoraCluster->GetIsolatedCaloHitList().end());

      pandora::FloatVector hitE, hitX, hitY, hitZ;
      edm4hep::MutableCluster     p_Cluster0 = pClusterCollection.create();
      edm4hep::MutableCluster*    p_Cluster  = &p_Cluster0;
      this->SetClusterSubDetectorEnergies(subDetectorNames, p_Cluster, pandoraCaloHitList, hitE, hitX, hitY, hitZ, calorimeterCollectionMap);

      float clusterCorrectEnergy(0.f);
      this->SetClusterEnergyAndError(pPandoraPfo, pPandoraCluster, p_Cluster, clusterCorrectEnergy);

      pandora::CartesianVector clusterPosition(0.f, 0.f, 0.f);
      const unsigned int       nHitsInCluster(pandoraCaloHitList.size());
      this->SetClusterPositionAndError(nHitsInCluster, hitE, hitX, hitY, hitZ, p_Cluster, clusterPosition);

      if (!hasTrack) {
        clustersWeightedPosition += clusterPosition * clusterCorrectEnergy;
        clustersTotalEnergy += clusterCorrectEnergy;
      }

      edm4hep::Cluster p_ClusterCon = *p_Cluster;
      pReconstructedParticle->addToClusters(p_ClusterCon);
    }

    if (!hasTrack) {
      if (clustersTotalEnergy < std::numeric_limits<float>::epsilon()) {
        MsgStream log(m_msgSvc, "PFOCreator");
        log << MSG::WARNING << "DDPfoCreator::CreateParticleFlowObjects: invalid cluster energy "
                               << clustersTotalEnergy << endmsg;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
      } else {
        referencePoint = clustersWeightedPosition * (1.f / clustersTotalEnergy);
      }
    } else {
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                               this->CalculateTrackBasedReferencePoint(pPandoraPfo, referencePoint, trackCollectionMap));
    }

    this->SetRecoParticleReferencePoint(referencePoint, pReconstructedParticle);
    this->AddTracksToRecoParticle(pPandoraPfo, pReconstructedParticle, trackCollectionMap);
    this->SetRecoParticlePropertiesFromPFO(pPandoraPfo, pReconstructedParticle);

    edm4hep::MutableVertex  pStartVertex0 = pStartVertexCollection.create();
    edm4hep::MutableVertex* pStartVertex  = &pStartVertex0;
    pStartVertex->setAlgorithmType(0);
    const float ref_value[3] = {referencePoint.GetX(), referencePoint.GetY(), referencePoint.GetZ()};
    pStartVertex->setPosition(edm4hep::Vector3f(ref_value));
    pStartVertex->addToParticles(*pReconstructedParticle);
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::InitialiseSubDetectorNames(pandora::StringVector& subDetectorNames) const {
  subDetectorNames.push_back("ecal");
  subDetectorNames.push_back("hcal");
  subDetectorNames.push_back("yoke");
  subDetectorNames.push_back("lcal");
  subDetectorNames.push_back("lhcal");
  subDetectorNames.push_back("bcal");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetClusterSubDetectorEnergies(const pandora::StringVector&   subDetectorNames,
                                                 edm4hep::MutableCluster* pCluster,
                                                 const pandora::CaloHitList&    pandoraCaloHitList,
                                                 pandora::FloatVector& hitE, pandora::FloatVector& hitX,
                                                 pandora::FloatVector& hitY, pandora::FloatVector& hitZ,
                                                 std::unordered_map<int, const edm4hep::CalorimeterHitCollection*> calorimeterCollections) const {
  
  std::vector<float> subDetectorEnergies;
  for (pandora::CaloHitList::const_iterator hIter = pandoraCaloHitList.begin(), hIterEnd = pandoraCaloHitList.end();
       hIter != hIterEnd; ++hIter) {
    const pandora::CaloHit* pPandoraCaloHit = *hIter;

    uint64_t storedID = reinterpret_cast<uint64_t>(pPandoraCaloHit->GetParentAddress());
    int collectionID = static_cast<int>(storedID >> 32);
    int index = static_cast<int>(storedID & 0xFFFFFFFF);
    auto it = calorimeterCollections.find(collectionID);
    if (it == calorimeterCollections.end()) {
        throw std::runtime_error("Collection ID not found!");
    }
    const edm4hep::CalorimeterHitCollection* collection = it->second;
    edm4hep::CalorimeterHit pCalorimeterHit = collection->at(index);

    pCluster->addToHits(pCalorimeterHit);

    const float caloHitEnergy(pCalorimeterHit.getEnergy());

    hitE.push_back(caloHitEnergy);
    hitX.push_back(pCalorimeterHit.getPosition().x);
    hitY.push_back(pCalorimeterHit.getPosition().y);
    hitZ.push_back(pCalorimeterHit.getPosition().z);

    subDetectorEnergies.resize(subDetectorNames.size());

    switch (CHT(pCalorimeterHit.getType()).caloID()) {
      case CHT::ecal:
        subDetectorEnergies[ECAL_INDEX] += caloHitEnergy;
        break;
      case CHT::hcal:
        subDetectorEnergies[HCAL_INDEX] += caloHitEnergy;
        break;
      case CHT::yoke:
        subDetectorEnergies[YOKE_INDEX] += caloHitEnergy;
        break;
      case CHT::lcal:
        subDetectorEnergies[LCAL_INDEX] += caloHitEnergy;
        break;
      case CHT::lhcal:
        subDetectorEnergies[LHCAL_INDEX] += caloHitEnergy;
        break;
      case CHT::bcal:
        subDetectorEnergies[BCAL_INDEX] += caloHitEnergy;
        break;
      default:
        MsgStream log(m_msgSvc, "PFOCreator");
        log << MSG::WARNING
            << "DDPfoCreator::SetClusterSubDetectorEnergies: no subdetector found for hit with type: "
            << pCalorimeterHit.getType() << endmsg;
    }
  }

  for (float subDetE : subDetectorEnergies) {
    pCluster->addToSubdetectorEnergies(subDetE);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetClusterEnergyAndError(const pandora::ParticleFlowObject* const pPandoraPfo,
                                            const pandora::Cluster* const            pPandoraCluster,
                                            edm4hep::MutableCluster* pCluster, float& clusterCorrectEnergy) const {
  const bool isEmShower((pandora::PHOTON == pPandoraPfo->GetParticleId()) ||
                        (pandora::E_MINUS == std::abs(pPandoraPfo->GetParticleId())));
  clusterCorrectEnergy = (isEmShower ? pPandoraCluster->GetCorrectedElectromagneticEnergy(m_pandora)
                                     : pPandoraCluster->GetCorrectedHadronicEnergy(m_pandora));

  if (clusterCorrectEnergy < std::numeric_limits<float>::epsilon())
    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

  const float stochasticTerm(isEmShower ? m_settings.m_emStochasticTerm : m_settings.m_hadStochasticTerm);
  const float constantTerm(isEmShower ? m_settings.m_emConstantTerm : m_settings.m_hadConstantTerm);
  const float energyError(
      std::sqrt(stochasticTerm * stochasticTerm / clusterCorrectEnergy + constantTerm * constantTerm) *
      clusterCorrectEnergy);

  pCluster->setEnergy(clusterCorrectEnergy);
  pCluster->setEnergyError(energyError);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetClusterPositionAndError(const unsigned int nHitsInCluster, pandora::FloatVector& hitE,
                                              pandora::FloatVector& hitX, pandora::FloatVector& hitY,
                                              pandora::FloatVector& hitZ, edm4hep::MutableCluster* pCluster,
                                              pandora::CartesianVector& clusterPositionVec) const {
  ClusterShapes* const pClusterShapes(
      new ClusterShapes(nHitsInCluster, hitE.data(), hitX.data(), hitY.data(), hitZ.data()));

  try {
    pCluster->setPhi(std::atan2(pClusterShapes->getEigenVecInertia()[1], pClusterShapes->getEigenVecInertia()[0]));
    pCluster->setITheta(std::acos(pClusterShapes->getEigenVecInertia()[2]));
    pCluster->setPosition(pClusterShapes->getCentreOfGravity());
    //ATTN these two lines below would only compile with ilcsoft HEAD V2015-10-13 and above
    //pLcioCluster->setPositionError(pClusterShapes->getCenterOfGravityErrors());
    //pLcioCluster->setDirectionError(pClusterShapes->getEigenVecInertiaErrors());
    clusterPositionVec.SetValues(pClusterShapes->getCentreOfGravity()[0], pClusterShapes->getCentreOfGravity()[1],
                                 pClusterShapes->getCentreOfGravity()[2]);
  } catch (...) {
    MsgStream log(m_msgSvc, "PFOCreator");
    log << MSG::WARNING << "DDPfoCreator::SetClusterPositionAndError: unidentified exception caught." << endmsg;
  }

  delete pClusterShapes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDPfoCreator::CalculateTrackBasedReferencePoint(
    const pandora::ParticleFlowObject* const pPandoraPfo, pandora::CartesianVector& referencePoint,
    std::unordered_map<int, const edm4hep::TrackCollection*> trackCollections) const {
  const pandora::TrackList& trackList(pPandoraPfo->GetTrackList());

  float                    totalTrackMomentumAtDca(0.f), totalTrackMomentumAtStart(0.f);
  pandora::CartesianVector referencePointAtDCAWeighted(0.f, 0.f, 0.f), referencePointAtStartWeighted(0.f, 0.f, 0.f);

  bool hasSiblings(false);
  for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end(); tIter != tIterEnd;
       ++tIter) {
    const pandora::Track* const pPandoraTrack(*tIter);

    if (!this->IsValidParentTrack(pPandoraTrack, trackList))
      continue;

    if (this->HasValidSiblingTrack(pPandoraTrack, trackList)) {
      // Presence of sibling tracks typically represents a conversion
      const pandora::CartesianVector& trackStartPoint((pPandoraTrack->GetTrackStateAtStart()).GetPosition());
      const float trackStartMomentum(((pPandoraTrack->GetTrackStateAtStart()).GetMomentum()).GetMagnitude());
      referencePointAtStartWeighted += trackStartPoint * trackStartMomentum;
      totalTrackMomentumAtStart += trackStartMomentum;
      hasSiblings = true;
    } else {
      uint64_t storedID = reinterpret_cast<uint64_t>(pPandoraTrack->GetParentAddress());
      int collectionID = static_cast<int>(storedID >> 32);
      int index = static_cast<int>(storedID & 0xFFFFFFFF);
      auto it = trackCollections.find(collectionID);
      if (it == trackCollections.end()) {
          throw std::runtime_error("Collection ID not found!");
      }
      const edm4hep::TrackCollection* collection = it->second;
      edm4hep::Track pTrack = collection->at(index);
    
      edm4hep::TrackState stateIP = pTrack.getTrackStates(edm4hep::TrackState::AtIP);

      const float              z0(pPandoraTrack->GetZ0());
      pandora::CartesianVector intersectionPoint(0.f, 0.f, 0.f);

      intersectionPoint.SetValues(stateIP.D0 * std::cos(stateIP.phi),
                                  stateIP.D0 * std::sin(stateIP.phi), z0);
      const float trackMomentumAtDca((pPandoraTrack->GetMomentumAtDca()).GetMagnitude());
      referencePointAtDCAWeighted += intersectionPoint * trackMomentumAtDca;
      totalTrackMomentumAtDca += trackMomentumAtDca;
    }
  }

  if (hasSiblings) {
    if (totalTrackMomentumAtStart < std::numeric_limits<float>::epsilon()) {
      MsgStream log(m_msgSvc, "PFOCreator");
      log << MSG::WARNING << "DDPfoCreator::CalculateTrackBasedReferencePoint: invalid track momentum "
                             << totalTrackMomentumAtStart << endmsg;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    } else {
      referencePoint = referencePointAtStartWeighted * (1.f / totalTrackMomentumAtStart);
    }
  } else {
    if (totalTrackMomentumAtDca < std::numeric_limits<float>::epsilon()) {
      MsgStream log(m_msgSvc, "PFOCreator");
      log << MSG::WARNING << "DDPfoCreator::CalculateTrackBasedReferencePoint: invalid track momentum "
                             << totalTrackMomentumAtDca << endmsg;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    } else {
      referencePoint = referencePointAtDCAWeighted * (1.f / totalTrackMomentumAtDca);
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::IsValidParentTrack(const pandora::Track* const pPandoraTrack,
                                      const pandora::TrackList&   allTrackList) const {
  const pandora::TrackList& parentTrackList(pPandoraTrack->GetParentList());

  for (pandora::TrackList::const_iterator iter = parentTrackList.begin(), iterEnd = parentTrackList.end();
       iter != iterEnd; ++iter) {
    if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
      continue;

    // ATTN This track must have a parent not in the all track list; still use it if it is the closest to the ip
    MsgStream log(m_msgSvc, "PFOCreator");
    log << MSG::WARNING
        << "DDPfoCreator::IsValidParentTrack: mismatch in track relationship information, use information as available "
        << endmsg;

    if (this->IsClosestTrackToIP(pPandoraTrack, allTrackList))
      return true;

    return false;
  }

  // Ideal case: All parents are associated to same pfo
  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::HasValidSiblingTrack(const pandora::Track* const pPandoraTrack,
                                        const pandora::TrackList&   allTrackList) const {
  const pandora::TrackList& siblingTrackList(pPandoraTrack->GetSiblingList());

  for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end();
       iter != iterEnd; ++iter) {
    if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
      continue;

    // ATTN This track must have a sibling not in the all track list; still use it if it has a second sibling that is in the list
    MsgStream log(m_msgSvc, "PFOCreator");
    log << MSG::WARNING << "DDPfoCreator::HasValidSiblingTrack: mismatch in track relationship information, use "
                              "information as available "
                           << endmsg;

    if (this->AreAnyOtherSiblingsInList(pPandoraTrack, allTrackList))
      return true;

    return false;
  }

  // Ideal case: All siblings associated to same pfo
  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::IsClosestTrackToIP(const pandora::Track* const pPandoraTrack,
                                      const pandora::TrackList&   allTrackList) const {
  const pandora::Track* pClosestTrack(NULL);
  float                 closestTrackDisplacement(std::numeric_limits<float>::max());

  for (pandora::TrackList::const_iterator iter = allTrackList.begin(), iterEnd = allTrackList.end(); iter != iterEnd;
       ++iter) {
    const pandora::Track* const pTrack(*iter);
    const float                 trialTrackDisplacement(pTrack->GetTrackStateAtStart().GetPosition().GetMagnitude());

    if (trialTrackDisplacement < closestTrackDisplacement) {
      closestTrackDisplacement = trialTrackDisplacement;
      pClosestTrack            = pTrack;
    }
  }

  return (pPandoraTrack == pClosestTrack);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDPfoCreator::AreAnyOtherSiblingsInList(const pandora::Track* const pPandoraTrack,
                                             const pandora::TrackList&   allTrackList) const {
  const pandora::TrackList& siblingTrackList(pPandoraTrack->GetSiblingList());

  for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end();
       iter != iterEnd; ++iter) {
    if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
      return true;
  }

  return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetRecoParticleReferencePoint(const pandora::CartesianVector&              referencePoint,
                                                 edm4hep::MutableReconstructedParticle* pReconstructedParticle) const {
  const float referencePointArray[3] = {referencePoint.GetX(), referencePoint.GetY(), referencePoint.GetZ()};
  pReconstructedParticle->setReferencePoint(referencePointArray);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::AddTracksToRecoParticle(const pandora::ParticleFlowObject* const     pPandoraPfo,
                                           edm4hep::MutableReconstructedParticle* pReconstructedParticle,
                                           std::unordered_map<int, const edm4hep::TrackCollection*> trackCollections) const {
  const pandora::TrackList& trackList(pPandoraPfo->GetTrackList());

  for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end(); tIter != tIterEnd;
       ++tIter) {
    const pandora::Track* const pTrack0(*tIter);
    uint64_t storedID = reinterpret_cast<uint64_t>(pTrack0->GetParentAddress());
    int collectionID = static_cast<int>(storedID >> 32);
    int index = static_cast<int>(storedID & 0xFFFFFFFF);
    auto it = trackCollections.find(collectionID);
    if (it == trackCollections.end()) {
        throw std::runtime_error("Collection ID not found!");
    }
    const edm4hep::TrackCollection* collection = it->second;
    edm4hep::Track pTrack = collection->at(index);

    MsgStream log(m_msgSvc, "PFOs");
    log << MSG::DEBUG << pTrack << endmsg;
    pReconstructedParticle->addToTracks(pTrack);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPfoCreator::SetRecoParticlePropertiesFromPFO(
    const pandora::ParticleFlowObject* const    pPandoraPfo,
    edm4hep::MutableReconstructedParticle *pReconstructedParticle) const {
  const float momentum[3] = {pPandoraPfo->GetMomentum().GetX(), pPandoraPfo->GetMomentum().GetY(),
                             pPandoraPfo->GetMomentum().GetZ()};
  pReconstructedParticle->setMomentum(momentum);
  pReconstructedParticle->setEnergy(pPandoraPfo->GetEnergy());
  pReconstructedParticle->setMass(pPandoraPfo->GetMass());
  pReconstructedParticle->setCharge(pPandoraPfo->GetCharge());
  pReconstructedParticle->setPDG(pPandoraPfo->GetParticleId());
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDPfoCreator::Settings::Settings()
    : m_emStochasticTerm(0.17f), m_hadStochasticTerm(0.6f), m_emConstantTerm(0.01f), m_hadConstantTerm(0.03f) {}
