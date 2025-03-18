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
 *  @file   DDMarlinPandora/src/DDMCParticleCreator.cc
 *
 *  @brief  Implementation of the mc particle creator class.
 *
 *  $Log: $
 */

#include "edm4hep/CaloHitContribution.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCRecoCaloAssociation.h"
#include "edm4hep/MCRecoTrackerAssociation.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackState.h"

#include "DDMCParticleCreator.h"

#include <string>

//forward declarations. See in DDPandoraPFANewProcessor.cc
double getFieldFromCompact();

DDMCParticleCreator::DDMCParticleCreator(const Settings& settings, const pandora::Pandora* const pPandora, IMessageSvc* msgSvc)
    : m_settings(settings), m_pandora(*pPandora), m_bField(getFieldFromCompact()), m_msgSvc(msgSvc) {}

//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::~DDMCParticleCreator() {}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateMCParticles(const std::vector<const edm4hep::MCParticleCollection*>& MCParticleCollections) {
  MsgStream log(m_msgSvc, "MCParticleCreator");
  for (int colIndex = 0; colIndex < MCParticleCollections.size(); colIndex++) {
    try {
      const edm4hep::MCParticleCollection* pMCParticleCollection = MCParticleCollections[colIndex];
      const int nElements(pMCParticleCollection->size());

      if (0 == nElements)
        continue;

      log << MSG::DEBUG << "Creating MCParticles particles" << endmsg;
      for (int i = 0, iMax = nElements; i < iMax; ++i) {
        try {
          edm4hep::MCParticle pMcParticle = pMCParticleCollection->at(i);

          PandoraApi::MCParticle::Parameters mcParticleParameters;
          mcParticleParameters.m_energy         = pMcParticle.getEnergy();
          mcParticleParameters.m_particleId     = pMcParticle.getPDG();
          mcParticleParameters.m_mcParticleType = pandora::MC_3D;

          
          uint64_t ID = pMcParticle.getObjectID().collectionID;
          ID = (ID << 32) | pMcParticle.getObjectID().index;
          mcParticleParameters.m_pParentAddress = reinterpret_cast<void*>(ID);     

          mcParticleParameters.m_momentum       = pandora::CartesianVector(
              pMcParticle.getMomentum().x, pMcParticle.getMomentum().y, pMcParticle.getMomentum().z);
          mcParticleParameters.m_vertex = pandora::CartesianVector(
              pMcParticle.getVertex().x, pMcParticle.getVertex().y, pMcParticle.getVertex().z);
          mcParticleParameters.m_endpoint = pandora::CartesianVector(
              pMcParticle.getEndpoint().x, pMcParticle.getEndpoint().y, pMcParticle.getEndpoint().z);

          PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                  PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

          // Create parent-daughter relationships
          for (edm4hep::MCParticle daughter : pMcParticle.getDaughters()) {
            uint64_t daughterID = daughter.getObjectID().collectionID;
            daughterID = (daughterID << 32) | daughter.getObjectID().index;    
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(
              m_pandora, reinterpret_cast<void*>(ID), reinterpret_cast<void*>(daughterID)));
          }
        } catch (pandora::StatusCodeException& statusCodeException) {
          log << MSG::ERROR << "Failed to extract MCParticle: " << statusCodeException.ToString() << endmsg;
        }
      }
    } catch(...) {
      log << MSG::ERROR << "Failed to extract MCParticle collection" << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateTrackToMCParticleRelationships(const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& linkCollections,
                                                                              const TrackVector&    trackVector,
                                                                              const std::vector<const edm4hep::TrackCollection*>& trackCollections) const {
  MsgStream log(m_msgSvc, "MCParticleCreator");
  std::unordered_map<int, const edm4hep::TrackCollection*> trackCollectionMap;
  for (auto* collection : trackCollections) { trackCollectionMap[collection->getID()] = collection; }
  
  for (unsigned ik = 0; ik < trackVector.size(); ik++) {
    ///////FINISH THIS////////
    uint64_t trackID = trackVector.at(ik);
    int collectionID = static_cast<int>(trackID >> 32);
    int index = static_cast<int>(trackID & 0xFFFFFFFF);
    auto it = trackCollectionMap.find(collectionID);
    if (it == trackCollectionMap.end()) {
        throw std::runtime_error("Collection ID not found!");
    }
    const edm4hep::TrackCollection* collection = it->second;
    edm4hep::Track pTrack = collection->at(index);

    // Get reconstructed momentum at dca
    edm4hep::TrackState stateAtIP = pTrack.getTrackStates(0);
    const pandora::Helix helixFit(stateAtIP.phi, stateAtIP.D0,
                                  stateAtIP.Z0, stateAtIP.omega,
                                  stateAtIP.tanLambda, m_bField);
    const float          recoMomentum(helixFit.GetMomentum().GetMagnitude());

    // Use momentum magnitude to identify best mc particle
    uint64_t pBestMCParticle = 0;
    float                bestDeltaMomentum(std::numeric_limits<float>::max());
    try {
      for (int colIndex = 0; colIndex < linkCollections.size(); colIndex++) {
        const edm4hep::TrackerHitSimTrackerHitLinkCollection *pLinkCollection = linkCollections[colIndex];
        for (unsigned ith = 0; ith < pTrack.trackerHits_size(); ith++) {
          for (unsigned ic = 0; ic < pLinkCollection->size(); ic++) {
            if ((pLinkCollection->at(ic)).getFrom() != pTrack.getTrackerHits(ith))
              continue;
            const edm4hep::SimTrackerHit pSimHit = (pLinkCollection->at(ic)).getTo();
            const edm4hep::MCParticle    ipa     = pSimHit.getParticle();

            const float trueMomentum(
                pandora::CartesianVector(ipa.getMomentum().x, ipa.getMomentum().y, ipa.getMomentum().z)
                    .GetMagnitude());
            const float deltaMomentum(std::fabs(recoMomentum - trueMomentum));
            if (deltaMomentum < bestDeltaMomentum) {
              uint64_t ID = ipa.getObjectID().collectionID;
              ID = (ID << 32) | ipa.getObjectID().index;
              pBestMCParticle   = ID;
              bestDeltaMomentum = deltaMomentum;
            }
          }
        }
      }

      if (0 == pBestMCParticle)
        continue;
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackToMCParticleRelationship(
        m_pandora, reinterpret_cast<void*>(trackID), reinterpret_cast<void*>(pBestMCParticle)));
    } catch (pandora::StatusCodeException& statusCodeException) {
      log << MSG::ERROR << "Failed to extract track to mc particle relationship: " << statusCodeException.ToString()
                           << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateCaloHitToMCParticleRelationships(const std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*>& linkCollections, 
                                                                                const CalorimeterHitVector& calorimeterHitVector,
                                                                                const std::vector<const edm4hep::CalorimeterHitCollection*>& eCalCollections,
                                                                                const std::vector<const edm4hep::CalorimeterHitCollection*>& hCalCollections,
                                                                                const std::vector<const edm4hep::CalorimeterHitCollection*>& mCalCollections,
                                                                                const std::vector<const edm4hep::CalorimeterHitCollection*>& lCalCollections,
                                                                                const std::vector<const edm4hep::CalorimeterHitCollection*>& lhCalCollections) const {
  MsgStream log(m_msgSvc, "MCParticleCreator");
  std::unordered_map<int, const edm4hep::CalorimeterHitCollection*> calorimeterCollectionMap;
  auto processCollections = [&](const std::vector<const edm4hep::CalorimeterHitCollection*>& collections) {
    for (auto* collection : collections) { calorimeterCollectionMap[collection->getID()] = collection; }
  };
  processCollections(eCalCollections);
  processCollections(hCalCollections);
  processCollections(mCalCollections);
  processCollections(lCalCollections);
  processCollections(lhCalCollections);

  typedef std::map<void*, float> MCParticleToEnergyWeightMap;
  MCParticleToEnergyWeightMap          mcParticleToEnergyWeightMap;

  for (int colIndex = 0; colIndex < linkCollections.size(); colIndex++) {
    try {
      const edm4hep::CaloHitSimCaloHitLinkCollection *pLinkCollection = linkCollections[colIndex];

      for (unsigned i_calo = 0; i_calo < calorimeterHitVector.size(); i_calo++) {
        try {
          uint64_t storedID = calorimeterHitVector.at(i_calo);
          int collectionID = static_cast<int>(storedID >> 32);
          int index = static_cast<int>(storedID & 0xFFFFFFFF);
          auto it = calorimeterCollectionMap.find(collectionID);
          if (it == calorimeterCollectionMap.end()) {
              throw std::runtime_error("Collection ID not found!");
          }
          const edm4hep::CalorimeterHitCollection* collection = it->second;
          edm4hep::CalorimeterHit pCalorimeterHit = collection->at(index);
          mcParticleToEnergyWeightMap.clear();

          for (unsigned ic = 0; ic < pLinkCollection->size(); ic++) {
            if ((pLinkCollection->at(ic)).getFrom() != pCalorimeterHit)
              continue;
            const edm4hep::SimCalorimeterHit pSimHit = (pLinkCollection->at(ic)).getTo();
            for (int iCont = 0, iEnd = pSimHit.contributions_size(); iCont < iEnd; ++iCont) {
              edm4hep::CaloHitContribution conb = pSimHit.getContributions(iCont);
              const edm4hep::MCParticle    ipa  = conb.getParticle();
              float                        ien  = conb.getEnergy();
              
              uint64_t ID =  ipa.getObjectID().collectionID;
              ID = (ID << 32) | ipa.getObjectID().index;
              void *mcVoidID = reinterpret_cast<void*>(ID);
              mcParticleToEnergyWeightMap[mcVoidID] += ien;
            }
          }

          for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter    = mcParticleToEnergyWeightMap.begin(),
                                                           mcParticleIterEnd = mcParticleToEnergyWeightMap.end();
               mcParticleIter != mcParticleIterEnd; ++mcParticleIter) {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(
              m_pandora, reinterpret_cast<void*>(storedID), mcParticleIter->first, mcParticleIter->second));
          }
        } catch (pandora::StatusCodeException& statusCodeException) {
          log << MSG::ERROR << "Failed to extract calo hit to mc particle relationship: "
                               << statusCodeException.ToString() << endmsg;
        } 
      }
    } catch(...) {
      log << MSG::ERROR << "Failed to extract Calo MCP Link collection" << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::Settings::Settings() {}
