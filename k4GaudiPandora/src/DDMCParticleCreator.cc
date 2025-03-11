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

#include "DDMCParticleCreator.h"

#include <string>

//forward declarations. See in DDPandoraPFANewProcessor.cc
double getFieldFromCompact();

DDMCParticleCreator::DDMCParticleCreator(const Settings& settings, const pandora::Pandora* const pPandora, IMessageSvc* msgSvc)
    : m_settings(settings), m_pandora(*pPandora), m_bField(getFieldFromCompact()), m_msgSvc(msgSvc) {
      m_id_pMC_map = std::make_unique<std::map<int, std::shared_ptr<edm4hep::MCParticle>>>();
    }

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
          std::shared_ptr<edm4hep::MCParticle> pMcParticle = std::make_shared<edm4hep::MCParticle>(pMCParticleCollection->at(i));
          m_id_pMC_map->emplace(pMcParticle->id().index, pMcParticle);
          //edm4hep::MCParticle pMcParticle0 = pMCParticleCollection->at(i);
	        //edm4hep::MCParticle *pMcParticle = &pMcParticle0;

          if (NULL == pMcParticle)
            log << MSG::ERROR << "Collection type mismatch" << endmsg;

          PandoraApi::MCParticle::Parameters mcParticleParameters;
          mcParticleParameters.m_energy         = pMcParticle->getEnergy();
          mcParticleParameters.m_particleId     = pMcParticle->getPDG();
          mcParticleParameters.m_mcParticleType = pandora::MC_3D;

          std::shared_ptr<void> pMCPVoidPtr = std::static_pointer_cast<void>(pMcParticle);
          void* pMCPVoid = pMCPVoidPtr.get();        
          mcParticleParameters.m_pParentAddress = pMCPVoid;

          mcParticleParameters.m_momentum       = pandora::CartesianVector(
              pMcParticle->getMomentum().x, pMcParticle->getMomentum().y, pMcParticle->getMomentum().z);
          mcParticleParameters.m_vertex = pandora::CartesianVector(
              pMcParticle->getVertex().x, pMcParticle->getVertex().y, pMcParticle->getVertex().z);
          mcParticleParameters.m_endpoint = pandora::CartesianVector(
              pMcParticle->getEndpoint().x, pMcParticle->getEndpoint().y, pMcParticle->getEndpoint().z);

          PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                  PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

          // Create parent-daughter relationships
          for (edm4hep::MCParticle daughter : pMcParticle->getDaughters()) {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                    PandoraApi::SetMCParentDaughterRelationship(m_pandora, pMcParticle.get(), &daughter));
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
                                                                              const TrackVector&    trackVector) const {
  MsgStream log(m_msgSvc, "MCParticleCreator");
  for (unsigned ik = 0; ik < trackVector.size(); ik++) {
    std::shared_ptr<edm4hep::Track> pTrack = trackVector.at(ik);
    // Get reconstructed momentum at dca
    const pandora::Helix helixFit(pTrack->getTrackStates(0).phi, pTrack->getTrackStates(0).D0,
                                  pTrack->getTrackStates(0).Z0, pTrack->getTrackStates(0).omega,
                                  pTrack->getTrackStates(0).tanLambda, m_bField);
    const float          recoMomentum(helixFit.GetMomentum().GetMagnitude());

    // Use momentum magnitude to identify best mc particle
    std::shared_ptr<edm4hep::MCParticle> pBestMCParticle = NULL;
    float                bestDeltaMomentum(std::numeric_limits<float>::max());
    try {
      for (int colIndex = 0; colIndex < linkCollections.size(); colIndex++) {
        const edm4hep::TrackerHitSimTrackerHitLinkCollection *pLinkCollection = linkCollections[colIndex];
        for (unsigned ith = 0; ith < pTrack->trackerHits_size(); ith++) {
          for (unsigned ic = 0; ic < pLinkCollection->size(); ic++) {
            if ((pLinkCollection->at(ic)).getFrom() != pTrack->getTrackerHits(ith))
              continue;
            const edm4hep::SimTrackerHit pSimHit = (pLinkCollection->at(ic)).getTo();
            const edm4hep::MCParticle    ipa     = pSimHit.getParticle();

            auto it = m_id_pMC_map->find(ipa.id().index);
            if (it == m_id_pMC_map->end())
              continue;
            const float trueMomentum(
                pandora::CartesianVector(ipa.getMomentum().x, ipa.getMomentum().y, ipa.getMomentum().z)
                    .GetMagnitude());
            const float deltaMomentum(std::fabs(recoMomentum - trueMomentum));
            if (deltaMomentum < bestDeltaMomentum) {
              pBestMCParticle   = it->second;
              bestDeltaMomentum = deltaMomentum;
            }
          }
        }
      }

      if (NULL == pBestMCParticle)
        continue;
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::SetTrackToMCParticleRelationship(m_pandora, pTrack.get(), pBestMCParticle.get()));
    } catch (pandora::StatusCodeException& statusCodeException) {
      log << MSG::ERROR << "Failed to extract track to mc particle relationship: " << statusCodeException.ToString()
                           << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateCaloHitToMCParticleRelationships(const std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*>& linkCollections, 
                                                                                const CalorimeterHitVector& calorimeterHitVector) const {
  MsgStream log(m_msgSvc, "MCParticleCreator");
  typedef std::map<const edm4hep::MCParticle*, float> MCParticleToEnergyWeightMap;
  MCParticleToEnergyWeightMap          mcParticleToEnergyWeightMap;

  for (int colIndex = 0; colIndex < linkCollections.size(); colIndex++) {
    try {
      const edm4hep::CaloHitSimCaloHitLinkCollection *pLinkCollection = linkCollections[colIndex];

      for (unsigned i_calo = 0; i_calo < calorimeterHitVector.size(); i_calo++) {
        try {
          mcParticleToEnergyWeightMap.clear();

          for (unsigned ic = 0; ic < pLinkCollection->size(); ic++) {
            if ((pLinkCollection->at(ic)).getFrom() != *(calorimeterHitVector.at(i_calo).get()))
              continue;
            const edm4hep::SimCalorimeterHit pSimHit = (pLinkCollection->at(ic)).getTo();
            for (int iCont = 0, iEnd = pSimHit.contributions_size(); iCont < iEnd; ++iCont) {
              edm4hep::CaloHitContribution conb = pSimHit.getContributions(iCont);
              const edm4hep::MCParticle    ipa  = conb.getParticle();
              float                        ien  = conb.getEnergy();
              
              auto it = m_id_pMC_map->find(ipa.id().index);
              if (it == m_id_pMC_map->end())
                continue;
              std::shared_ptr<edm4hep::MCParticle> p_tmp = it->second;
              mcParticleToEnergyWeightMap[p_tmp.get()] += ien;
            }
          }

          for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter    = mcParticleToEnergyWeightMap.begin(),
                                                           mcParticleIterEnd = mcParticleToEnergyWeightMap.end();
               mcParticleIter != mcParticleIterEnd; ++mcParticleIter) {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                    PandoraApi::SetCaloHitToMCParticleRelationship(
                                        m_pandora, calorimeterHitVector.at(i_calo).get(), mcParticleIter->first, mcParticleIter->second));
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
