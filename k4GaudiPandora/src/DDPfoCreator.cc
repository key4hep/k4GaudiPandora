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

#include "DDPfoCreator.h"

#include "CalorimeterHitType.h"
#include "ClusterShapes.h"

#include "Api/PandoraApi.h"

#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/Cluster.h>
#include <edm4hep/ClusterCollection.h>
#include <edm4hep/MutableCluster.h>
#include <edm4hep/ReconstructedParticle.h>
#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/Vector3f.h>
#include <edm4hep/Vertex.h>
#include <edm4hep/VertexCollection.h>

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>

DDPfoCreator::DDPfoCreator(const Settings& settings, pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm)
    : m_settings(settings), m_pandora(pandora), m_algorithm(*algorithm) {}

pandora::StatusCode
DDPfoCreator::CreateParticleFlowObjects(edm4hep::ClusterCollection& pClusterCollection,
                                        edm4hep::ReconstructedParticleCollection& pReconstructedParticleCollection,
                                        edm4hep::VertexCollection& pStartVertexCollection) const {
  const pandora::PfoList* pPandoraPfoList = NULL;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(m_pandora, pPandoraPfoList))

  pandora::StringVector subDetectorNames;
  this->InitialiseSubDetectorNames(subDetectorNames);
  // TODO: Add metadata
  // pClusterCollection->parameters().setValues("ClusterSubdetectorNames", subDetectorNames);

  // Create lcio "reconstructed particles" from the pandora "particle flow objects"
  for (const pandora::ParticleFlowObject* const pPandoraPfo : *pPandoraPfoList) {
    auto reconstructedParticle = pReconstructedParticleCollection.create();

    const bool hasTrack = !pPandoraPfo->GetTrackList().empty();

    float clustersTotalEnergy = 0.f;
    pandora::CartesianVector referencePoint{0.f, 0.f, 0.f}, clustersWeightedPosition{0.f, 0.f, 0.f};

    const pandora::ClusterList& clusterList(pPandoraPfo->GetClusterList());
    for (const auto* pPandoraCluster : clusterList) {
      pandora::CaloHitList pandoraCaloHitList;
      pPandoraCluster->GetOrderedCaloHitList().FillCaloHitList(pandoraCaloHitList);
      pandoraCaloHitList.insert(pandoraCaloHitList.end(), pPandoraCluster->GetIsolatedCaloHitList().begin(),
                                pPandoraCluster->GetIsolatedCaloHitList().end());

      pandora::FloatVector hitE, hitX, hitY, hitZ;
      auto cluster = pClusterCollection.create();
      setClusterSubDetectorEnergies(subDetectorNames, cluster, pandoraCaloHitList, hitE, hitX, hitY, hitZ);

      float clusterCorrectEnergy = 0.f;
      setClusterEnergyAndError(pPandoraPfo, pPandoraCluster, cluster, clusterCorrectEnergy);

      pandora::CartesianVector clusterPosition{0.f, 0.f, 0.f};
      setClusterPositionAndError(pandoraCaloHitList.size(), hitE, hitX, hitY, hitZ, cluster, clusterPosition);

      if (!hasTrack) {
        clustersWeightedPosition += clusterPosition * clusterCorrectEnergy;
        clustersTotalEnergy += clusterCorrectEnergy;
      }

      reconstructedParticle.addToClusters(cluster);
    }

    if (!hasTrack) {
      if (clustersTotalEnergy < std::numeric_limits<float>::epsilon()) {
        m_algorithm.warning() << "DDPfoCreator::CreateParticleFlowObjects: invalid cluster energy "
                              << clustersTotalEnergy << endmsg;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
      } else {
        referencePoint = clustersWeightedPosition * (1.f / clustersTotalEnergy);
      }
    } else {
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                               this->calculateTrackBasedReferencePoint(pPandoraPfo, referencePoint));
    }

    reconstructedParticle.setReferencePoint({referencePoint.GetX(), referencePoint.GetY(), referencePoint.GetZ()});
    this->AddTracksToRecoParticle(pPandoraPfo, reconstructedParticle);
    this->SetRecoParticlePropertiesFromPFO(pPandoraPfo, reconstructedParticle);

    auto startVertex = pStartVertexCollection.create();
    startVertex.setAlgorithmType(0);
    startVertex.setPosition({referencePoint.GetX(), referencePoint.GetY(), referencePoint.GetZ()});
    startVertex.addToParticles(reconstructedParticle);

    reconstructedParticle.setDecayVertex(startVertex);
  }

  return pandora::STATUS_CODE_SUCCESS;
}

void DDPfoCreator::InitialiseSubDetectorNames(pandora::StringVector& subDetectorNames) const {
  subDetectorNames.push_back("ecal");
  subDetectorNames.push_back("hcal");
  subDetectorNames.push_back("yoke");
  subDetectorNames.push_back("lcal");
  subDetectorNames.push_back("lhcal");
  subDetectorNames.push_back("bcal");
}

void DDPfoCreator::setClusterSubDetectorEnergies(const pandora::StringVector& subDetectorNames,
                                                 edm4hep::MutableCluster& cluster,
                                                 const pandora::CaloHitList& pandoraCaloHitList,
                                                 pandora::FloatVector& hitE, pandora::FloatVector& hitX,
                                                 pandora::FloatVector& hitY, pandora::FloatVector& hitZ) const {
  std::vector<float> subDetectorEnergies;
  subDetectorEnergies.reserve(subDetectorNames.size());
  for (size_t i = 0; i < subDetectorNames.size() && i < cluster.getSubdetectorEnergies().size(); ++i) {
    subDetectorEnergies.push_back(cluster.getSubdetectorEnergies()[i]);
  }
  for (size_t i = cluster.getSubdetectorEnergies().size(); i < subDetectorNames.size(); ++i) {
    subDetectorEnergies.push_back(0.f);
  }
  for (const auto* pPandoraCaloHit : pandoraCaloHitList) {
    const auto& pCalorimeterHit = *static_cast<const edm4hep::CalorimeterHit*>(pPandoraCaloHit->GetParentAddress());

    cluster.addToHits(pCalorimeterHit);

    const float caloHitEnergy = pCalorimeterHit.getEnergy();
    hitE.push_back(caloHitEnergy);
    hitX.push_back(pCalorimeterHit.getPosition()[0]);
    hitY.push_back(pCalorimeterHit.getPosition()[1]);
    hitZ.push_back(pCalorimeterHit.getPosition()[2]);

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
      m_algorithm.warning()
          << "DDPfoCreator::setClusterSubDetectorEnergies: no subdetector found for hit with type: "
          << pCalorimeterHit.getType() << endmsg;
      break; // pass
    }
  }
  for (size_t i = 0; i < subDetectorEnergies.size(); ++i) {
    cluster.addToSubdetectorEnergies(subDetectorEnergies[i]);
  }
}

void DDPfoCreator::setClusterEnergyAndError(const pandora::ParticleFlowObject* const pPandoraPfo,
                                            const pandora::Cluster* const pPandoraCluster,
                                            edm4hep::MutableCluster& cluster, float& clusterCorrectEnergy) const {
  const bool isEmShower((pandora::PHOTON == pPandoraPfo->GetParticleId()) ||
                        (pandora::E_MINUS == std::abs(pPandoraPfo->GetParticleId())));
  clusterCorrectEnergy = (isEmShower ? pPandoraCluster->GetCorrectedElectromagneticEnergy(m_pandora)
                                     : pPandoraCluster->GetCorrectedHadronicEnergy(m_pandora));

  if (clusterCorrectEnergy < std::numeric_limits<float>::epsilon())
    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

  const float stochasticTerm = isEmShower ? m_settings.m_emStochasticTerm : m_settings.m_hadStochasticTerm;
  const float constantTerm = isEmShower ? m_settings.m_emConstantTerm : m_settings.m_hadConstantTerm;
  const float energyError =
      std::sqrt(stochasticTerm * stochasticTerm / clusterCorrectEnergy + constantTerm * constantTerm) *
      clusterCorrectEnergy;

  cluster.setEnergy(clusterCorrectEnergy);
  cluster.setEnergyError(energyError);
}

void DDPfoCreator::setClusterPositionAndError(const std::size_t nHitsInCluster, pandora::FloatVector& hitE,
                                              pandora::FloatVector& hitX, pandora::FloatVector& hitY,
                                              pandora::FloatVector& hitZ, edm4hep::MutableCluster& cluster,
                                              pandora::CartesianVector& clusterPositionVec) const {
  ClusterShapes clusterShape(nHitsInCluster, hitE.data(), hitX.data(), hitY.data(), hitZ.data());

  try {
    cluster.setPhi(std::atan2(clusterShape.getEigenVecInertia()[1], clusterShape.getEigenVecInertia()[0]));
    cluster.setITheta(std::acos(clusterShape.getEigenVecInertia()[2]));
    cluster.setPosition(clusterShape.getCentreOfGravity());
    // TODO: Enable these that are not originally in DDMarlinPandora
    // cluster.setPositionError(std::array<float, 6>{clusterShape.getCenterOfGravityErrors()});
    // cluster.setDirectionError(clusterShape.getEigenVecInertiaErrors());
    clusterPositionVec.SetValues(clusterShape.getCentreOfGravity()[0], clusterShape.getCentreOfGravity()[1],
                                 clusterShape.getCentreOfGravity()[2]);
  } catch (...) {
    m_algorithm.warning() << "DDPfoCreator::setClusterPositionAndError: unidentified exception caught." << endmsg;
  }
}

pandora::StatusCode
DDPfoCreator::calculateTrackBasedReferencePoint(const pandora::ParticleFlowObject* const pPandoraPfo,
                                                pandora::CartesianVector& referencePoint) const {
  const pandora::TrackList& trackList = pPandoraPfo->GetTrackList();

  float totalTrackMomentumAtDca = 0.f, totalTrackMomentumAtStart = 0.f;
  pandora::CartesianVector referencePointAtDCAWeighted{0.f, 0.f, 0.f}, referencePointAtStartWeighted{0.f, 0.f, 0.f};

  bool hasSiblings = false;
  for (const pandora::Track* pPandoraTrack : trackList) {
    if (!this->isValidParentTrack(pPandoraTrack, trackList))
      continue;

    if (this->hasValidSiblingTrack(pPandoraTrack, trackList)) {
      // Presence of sibling tracks typically represents a conversion
      const pandora::CartesianVector& trackStartPoint((pPandoraTrack->GetTrackStateAtStart()).GetPosition());
      const float trackStartMomentum(((pPandoraTrack->GetTrackStateAtStart()).GetMomentum()).GetMagnitude());
      referencePointAtStartWeighted += trackStartPoint * trackStartMomentum;
      totalTrackMomentumAtStart += trackStartMomentum;
      hasSiblings = true;
    } else {
      const float z0(pPandoraTrack->GetZ0());
      pandora::CartesianVector intersectionPoint(0.f, 0.f, 0.f);

      const auto& pTrack = *static_cast<const edm4hep::Track*>(pPandoraTrack->GetParentAddress());
      auto const& trackState = pTrack.getTrackStates(0);

      intersectionPoint.SetValues(trackState.D0 * std::cos(trackState.phi), trackState.D0 * std::sin(trackState.phi),
                                  z0);
      const float trackMomentumAtDca = pPandoraTrack->GetMomentumAtDca().GetMagnitude();
      referencePointAtDCAWeighted += intersectionPoint * trackMomentumAtDca;
      totalTrackMomentumAtDca += trackMomentumAtDca;
    }
  }

  if (hasSiblings) {
    if (totalTrackMomentumAtStart < std::numeric_limits<float>::epsilon()) {
      m_algorithm.warning() << "DDPfoCreator::calculateTrackBasedReferencePoint: invalid track momentum "
                            << totalTrackMomentumAtStart << endmsg;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    } else {
      referencePoint = referencePointAtStartWeighted * (1.f / totalTrackMomentumAtStart);
    }
  } else {
    if (totalTrackMomentumAtDca < std::numeric_limits<float>::epsilon()) {
      m_algorithm.warning() << "DDPfoCreator::calculateTrackBasedReferencePoint: invalid track momentum "
                            << totalTrackMomentumAtDca << endmsg;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    } else {
      referencePoint = referencePointAtDCAWeighted * (1.f / totalTrackMomentumAtDca);
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

bool DDPfoCreator::isValidParentTrack(const pandora::Track* const pPandoraTrack,
                                      const pandora::TrackList& allTrackList) const {
  const pandora::TrackList& parentTrackList(pPandoraTrack->GetParentList());

  for (const auto& track : parentTrackList) {
    if (std::ranges::find(allTrackList, track) != allTrackList.end())
      continue;

    // ATTN This track must have a parent not in the all track list; still use it if it is the closest to the ip
    m_algorithm.warning()
        << "DDPfoCreator::isValidParentTrack: mismatch in track relationship information, use information as available "
        << endmsg;

    if (isClosestTrackToIP(pPandoraTrack, allTrackList))
      return true;

    return false;
  }

  // Ideal case: All parents are associated to same pfo
  return true;
}

bool DDPfoCreator::hasValidSiblingTrack(const pandora::Track* const pPandoraTrack,
                                        const pandora::TrackList& allTrackList) const {
  const pandora::TrackList& siblingTrackList(pPandoraTrack->GetSiblingList());

  for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end();
       iter != iterEnd; ++iter) {
    if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
      continue;

    // ATTN This track must have a sibling not in the all track list; still use it if it has a second sibling that is in
    // the list
    m_algorithm.warning()
        << "DDPfoCreator::hasValidSiblingTrack: mismatch in track relationship information, use information as available "
        << endmsg;

    if (this->AreAnyOtherSiblingsInList(pPandoraTrack, allTrackList))
      return true;

    return false;
  }

  // Ideal case: All siblings associated to same pfo
  return true;
}

bool DDPfoCreator::isClosestTrackToIP(const pandora::Track* const pPandoraTrack,
                                      const pandora::TrackList& allTrackList) const {
  const pandora::Track* pClosestTrack = nullptr;
  float closestTrackDisplacement = std::numeric_limits<float>::max();

  for (const pandora::Track* const pTrack : allTrackList) {
    const float trialTrackDisplacement(pTrack->GetTrackStateAtStart().GetPosition().GetMagnitude());

    if (trialTrackDisplacement < closestTrackDisplacement) {
      closestTrackDisplacement = trialTrackDisplacement;
      pClosestTrack = pTrack;
    }
  }

  return pPandoraTrack == pClosestTrack;
}

bool DDPfoCreator::AreAnyOtherSiblingsInList(const pandora::Track* const pPandoraTrack,
                                             const pandora::TrackList& allTrackList) const {
  for (const auto& track : pPandoraTrack->GetSiblingList()) {
    if (std::ranges::find(allTrackList, track) != allTrackList.end())
      return true;
  }

  return false;
}

void DDPfoCreator::AddTracksToRecoParticle(const pandora::ParticleFlowObject* const pPandoraPfo,
                                           edm4hep::MutableReconstructedParticle& pReconstructedParticle) const {
  for (const auto* pTrack : pPandoraPfo->GetTrackList()) {
    const auto& pLcioTrack = *static_cast<const edm4hep::Track*>(pTrack->GetParentAddress());
    pReconstructedParticle.addToTracks(pLcioTrack);
  }
}

void DDPfoCreator::SetRecoParticlePropertiesFromPFO(
    const pandora::ParticleFlowObject* const pPandoraPfo,
    edm4hep::MutableReconstructedParticle& reconstructedParticle) const {
  const float momentum[3] = {pPandoraPfo->GetMomentum().GetX(), pPandoraPfo->GetMomentum().GetY(),
                             pPandoraPfo->GetMomentum().GetZ()};
  reconstructedParticle.setMomentum(momentum);
  reconstructedParticle.setEnergy(pPandoraPfo->GetEnergy());
  reconstructedParticle.setMass(pPandoraPfo->GetMass());
  reconstructedParticle.setCharge(pPandoraPfo->GetCharge());
  reconstructedParticle.setPDG(pPandoraPfo->GetParticleId());
}

DDPfoCreator::Settings::Settings()
    : m_emStochasticTerm(0.17f), m_hadStochasticTerm(0.6f), m_emConstantTerm(0.01f), m_hadConstantTerm(0.03f) {}

