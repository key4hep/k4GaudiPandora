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
#include "PandoraPFAlg.h"

#include <string>

//forward declarations. See in DDPandoraPFANewProcessor.cc
double getFieldFromCompact();

DDMCParticleCreator::DDMCParticleCreator(const Settings& settings, const pandora::Pandora* const pPandora, MsgStream log)
    : m_settings(settings), m_pandora(*pPandora), m_bField(getFieldFromCompact()), m_log(log) {}

//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::~DDMCParticleCreator() {}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateMCParticles(std::vector<const edm4hep::MCParticleCollection*>& MCParticleCollections) const {
  for (int colIndex = 0; colIndex < MCParticleCollections.size(); colIndex++) {
    try {
      const edm4hep::CalorimeterHitCollection* pMCParticleCollection = MCParticleCollections[colIndex];
      const int nElements(pMCParticleCollection->size());

      if (0 == nElements)
        continue;

      m_log << MSG::DEBUG << "Creating " << m_settings.m_mcParticleCollections[colIndex] << " particles" << endmsg;
      for (int i = 0, iMax = pMCParticleCollection->size(); i < iMax; ++i) {
        try {
          edm4hep::MCParticle pMcParticle = dynamic_cast<MCParticle>(pMCParticleCollection->at(i));

          if (NULL == pMcParticle)
            m_log << MSG::ERROR << "Collection type mismatch" << endmsg;

          PandoraApi::MCParticle::Parameters mcParticleParameters;
          mcParticleParameters.m_energy         = pMcParticle.getEnergy();
          mcParticleParameters.m_particleId     = pMcParticle.getPDG();
          mcParticleParameters.m_mcParticleType = pandora::MC_3D;
          mcParticleParameters.m_pParentAddress = pMcParticle;
          mcParticleParameters.m_momentum       = pandora::CartesianVector(
              pMcParticle.getMomentum().x, pMcParticle.getMomentum().y, pMcParticle.getMomentum().z);
          mcParticleParameters.m_vertex = pandora::CartesianVector(
              pMcParticle.getVertex().x, pMcParticle.getVertex().y, pMcParticle.getVertex().z);
          mcParticleParameters.m_endpoint = pandora::CartesianVector(
              pMcParticle.getEndpoint().x, pMcParticle.getEndpoint().y, pMcParticle.getEndpoint().z);

          PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                  PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

          // Create parent-daughter relationships
          for (MCParticleVec::const_iterator itDaughter    = pMcParticle.getDaughters().begin(),
                                             itDaughterEnd = pMcParticle.getDaughters().end();
               itDaughter != itDaughterEnd; ++itDaughter) {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                    PandoraApi::SetMCParentDaughterRelationship(m_pandora, pMcParticle, *itDaughter));
          }
        } catch (pandora::StatusCodeException& statusCodeException) {
          m_log << MSG::ERROR << "Failed to extract MCParticle: " << statusCodeException.ToString() << endmsg;
        }
      }
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateTrackToMCParticleRelationships(std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& linkCollections,
                                                                              const TrackVector&    trackVector) const {
  for (unsigned ik = 0; ik < trackVector.size(); ik++) {
    const edm4hep::Track* pTrack = trackVector.at(ik);
    // Get reconstructed momentum at dca
    const pandora::Helix helixFit(pTrack->getTrackStates(0).phi, pTrack->getTrackStates(0).D0,
                                  pTrack->getTrackStates(0).Z0, pTrack->getTrackStates(0).omega,
                                  pTrack->getTrackStates(0).tanLambda, m_bField);
    const float          recoMomentum(helixFit.GetMomentum().GetMagnitude());

    // Use momentum magnitude to identify best mc particle
    edm4hep::MCParticle* pBestMCParticle = NULL;
    float                bestDeltaMomentum(std::numeric_limits<float>::max());
    try {
      for (int colIndex = 0; colIndex < linkCollections.size(); colIndex++) {
        const edm4hep::TrackerHitSimTrackerHitLinkCollection pLinkCollection = linkCollections[colIndex];
        for (unsigned ith = 0; ith < pTrack->trackerHits_size(); ith++) {
          for (unsigned ic = 0; ic < pLinkCollection.size(); ic++) {
            if (pLinkCollection.at(ic).getFrom() != pTrack->getTrackerHits(ith))
              continue;
            const edm4hep::SimTrackerHit pSimHit = pLinkCollection.at(ic).getTo();
            const edm4hep::MCParticle    ipa     = pSimHit.getParticle();
            if (m_id_pMC_map->find(ipa.id()) == m_id_pMC_map->end())
              continue;
            const float trueMomentum(
                pandora::CartesianVector(ipa.getMomentum().x, ipa.getMomentum().y, ipa.getMomentum().z)
                    .GetMagnitude());
            const float deltaMomentum(std::fabs(recoMomentum - trueMomentum));
            if (deltaMomentum < bestDeltaMomentum) {
              pBestMCParticle   = const_cast<edm4hep::MCParticle*>((*m_id_pMC_map)[ipa.id()]);
              bestDeltaMomentum = deltaMomentum;
            }
          }
        }
      }

      if (NULL == pBestMCParticle)
        continue;
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::SetTrackToMCParticleRelationship(m_pandora, pTrack, pBestMCParticle));
    } catch (pandora::StatusCodeException& statusCodeException) {
      m_log << MSG::ERROR << "Failed to extract track to mc particle relationship: " << statusCodeException.ToString()
                           << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateCaloHitToMCParticleRelationships(std::vector<const edm4hep::CaloHitMCParticleLinkCollection*>& linkCollections, 
                                                                                const CalorimeterHitVector& calorimeterHitVector) const {
  typedef std::map<MCParticle*, float> MCParticleToEnergyWeightMap;
  MCParticleToEnergyWeightMap          mcParticleToEnergyWeightMap;

  for (int colIndex = 0; colIndex < linkCollections.size(); colIndex++) {
    try {
      const std::vector<edm4hep::CaloHitMCParticleLinkCollection>& pLinkCollection = linkCollections[colIndex];

      for (unsigned i_calo = 0; i_calo < calorimeterHitVector.size(); i_calo++) {
        try {
          mcParticleToEnergyWeightMap.clear();

          for (unsigned ic = 0; ic < pLinkCollection.size(); ic++) {
            if (pLinkCollection.at(ic).getFrom() != (*(calorimeterHitVector.at(i_calo))))
              continue;
            const edm4hep::SimCalorimeterHit pSimHit = pLinkCollection.at(ic).getTo();
            for (int iCont = 0, iEnd = pSimHit.contributions_size(); iCont < iEnd; ++iCont) {
              edm4hep::CaloHitContribution conb = pSimHit.getContributions(iCont);
              const edm4hep::MCParticle    ipa  = conb.getParticle();
              float                        ien  = conb.getEnergy();
              if (m_id_pMC_map->find(ipa.id()) == m_id_pMC_map->end())
                continue;
              const edm4hep::MCParticle* p_tmp = (*m_id_pMC_map)[ipa.id()];
              mcParticleToEnergyWeightMap[p_tmp] += ien;
            }
          }

          for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter    = mcParticleToEnergyWeightMap.begin(),
                                                           mcParticleIterEnd = mcParticleToEnergyWeightMap.end();
               mcParticleIter != mcParticleIterEnd; ++mcParticleIter) {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                    PandoraApi::SetCaloHitToMCParticleRelationship(
                                        m_pandora, *(calorimeterHitVector.at(i_calo)), mcParticleIter->first, mcParticleIter->second));
          }
        } catch (pandora::StatusCodeException& statusCodeException) {
          m_log << MSG::ERROR << "Failed to extract calo hit to mc particle relationship: "
                               << statusCodeException.ToString() << endmsg;
        } 
      }
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::Settings::Settings() {}
