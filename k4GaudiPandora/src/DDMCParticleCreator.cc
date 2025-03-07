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

#include "Gaudi/Property.h"
#include "GaudiAlg/Transformer.h"
// Define BaseClass_t
#include "edm4hep/CaloHitContribution.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCRecoCaloAssociation.h"
#include "edm4hep/MCRecoTrackerAssociation.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/Track.h"
#include "k4FWCore/BaseClass.h"

#include "DDMCParticleCreator.h"
#include "PandoraPFAlg.h"

#include <string>

//forward declarations. See in DDPandoraPFANewProcessor.cc
double getFieldFromCompact();

DDMCParticleCreator::DDMCParticleCreator(const Settings& settings, const pandora::Pandora* const pPandora)
    : m_settings(settings), m_pandora(*pPandora), m_bField(getFieldFromCompact()) {}

//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::~DDMCParticleCreator() {}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateMCParticles(
    std::vector<const edm4hep::MCParticleCollection&> mcParticleCollections) const {
  try {
    for (const auto& mcParticleCollection : mcParticleCollections) {
      if (mcParticleCollection.empty()) {
        throw std::runtime_error("MCParticle collection is empty. No particles to process.");
      }

      for (const auto& mcParticle : mcParticleCollection) {
        try {
          PandoraApi::MCParticle::Parameters mcParticleParameters;
          mcParticleParameters.m_energy         = mcParticle.getEnergy();
          mcParticleParameters.m_particleId     = mcParticle.getPDG();
          mcParticleParameters.m_mcParticleType = pandora::MC_3D;
          mcParticleParameters.m_pParentAddress = &mcParticle;
          mcParticleParameters.m_momentum = pandora::CartesianVector(
          mcParticle.getMomentum().x, mcParticle.getMomentum().y, mcParticle.getMomentum().z);
          mcParticleParameters.m_vertex = pandora::CartesianVector(
          mcParticle.getVertex().x, mcParticle.getVertex().y, mcParticle.getVertex().z);
          mcParticleParameters.m_endpoint = pandora::CartesianVector(
          mcParticle.getEndpoint().x, mcParticle.getEndpoint().y, mcParticle.getEndpoint().z);

          PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                  PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

          // Create parent-daughter relationships
          for (const auto& daughter : mcParticle.getDaughters()) {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                    PandoraApi::SetMCParentDaughterRelationship(m_pandora, &mcParticle, &daughter));
          }
        } catch (const pandora::StatusCodeException& statusCodeException) {
          streamlog_out(ERROR) << "Failed to extract MCParticle: " << statusCodeException.ToString() << std::endl;
        } catch (const std::exception& exception) {
          streamlog_out(WARNING) << "Failed to extract MCParticle: " << exception.what() << std::endl;
        }
      }
    }
  } catch (const std::exception& exception) {
    streamlog_out(DEBUG5) << "Failed to extract MCParticles collections: " << exception.what() << std::endl;
    throw;
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateTrackToMCParticleRelationships(
    const CollectionMaps& collectionMaps, const TrackVector& trackVector) const {
  
  for (const auto& pTrack : trackVector) {
    // Get reconstructed momentum at dca
    const auto& trackState = pTrack->getTrackStates(0);
    const pandora::Helix helixFit(trackState.phi, trackState.D0, trackState.Z0, trackState.omega,
                                  trackState.tanLambda, m_bField);
    const float recoMomentum = helixFit.GetMomentum().GetMagnitude();

    // Use momentum magnitude to identify best MC particle
    edm4hep::MCParticle* pBestMCParticle = nullptr;
    float bestDeltaMomentum = std::numeric_limits<float>::max();

    try {
      for (const auto& collectionName : m_settings.m_TrackRelationCollections) {
        auto it = collectionMaps.collectionMap_TrkRel.find(collectionName);
        if (it == collectionMaps.collectionMap_TrkRel.end()) {
          continue;
        }

        // Retrieve MCRecoTrackerAssociation collection
        const std::vector<edm4hep::MCRecoTrackerAssociation>& pMCRecoTrackerAssociationCollection = it->second;

        for (unsigned ith = 0; ith < pTrack->trackerHits_size(); ++ith) {
          for (unsigned ic = 0; ic < pMCRecoTrackerAssociationCollection.size(); ++ic) {
            if (pMCRecoTrackerAssociationCollection.at(ic).getRec().id() != pTrack->getTrackerHits(ith).id())
              continue;

            const edm4hep::ConstSimTrackerHit pSimHit = pMCRecoTrackerAssociationCollection.at(ic).getSim();
            const edm4hep::ConstMCParticle ipa = pSimHit.getMCParticle();

            if (m_id_pMC_map->find(ipa.id()) == m_id_pMC_map->end())
              continue;

            const float trueMomentum =
                pandora::CartesianVector(ipa.getMomentum()[0], ipa.getMomentum()[1], ipa.getMomentum()[2])
                    .GetMagnitude();
            const float deltaMomentum = std::fabs(recoMomentum - trueMomentum);

            if (deltaMomentum < bestDeltaMomentum) {
              pBestMCParticle = const_cast<edm4hep::MCParticle*>((*m_id_pMC_map)[ipa.id()]);
              bestDeltaMomentum = deltaMomentum;
            }
          }
        }
      }

      if (pBestMCParticle == nullptr) {
        streamlog_out(WARNING) << "No suitable MC particle found for track association." << std::endl;
        continue;
      }

      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::SetTrackToMCParticleRelationship(m_pandora, pTrack, pBestMCParticle));

    } catch (const pandora::StatusCodeException& statusCodeException) {
      streamlog_out(ERROR) << "Failed to extract track to MC particle relationship: "
                           << statusCodeException.ToString() << std::endl;
    } catch (const std::exception& exception) {
      streamlog_out(WARNING) << "Exception encountered: " << exception.what() << std::endl;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateCaloHitToMCParticleRelationships(
    const CollectionMaps& collectionMaps, const std::vector<edm4hep::CalorimeterHit>& calorimeterHitVector) const {
  using MCParticleToEnergyWeightMap = std::map<const edm4hep::MCParticle*, float>;
  MCParticleToEnergyWeightMap mcParticleToEnergyWeightMap;

  for (StringVector::const_iterator iter = m_settings.m_lcCaloHitRelationCollections.begin(),
                                    iterEnd = m_settings.m_lcCaloHitRelationCollections.end();
       iter != iterEnd; ++iter) {
    if (collectionMaps.collectionMap_CaloRel.find(*iter) == collectionMaps.collectionMap_CaloRel.end())
      continue;

    try {
      const std::vector<edm4hep::MCRecoCaloAssociation>& pMCRecoCaloAssociationCollection =
          (collectionMaps.collectionMap_CaloRel.find(*iter))->second;

      for (size_t i_calo = 0; i_calo < calorimeterHitVector.size(); i_calo++) {
        try {
          mcParticleToEnergyWeightMap.clear();

          for (size_t ic = 0; ic < pMCRecoCaloAssociationCollection.size(); ic++) {
            if (pMCRecoCaloAssociationCollection.at(ic).getRec().id() != calorimeterHitVector.at(i_calo).id())
              continue;

            const edm4hep::ConstSimCalorimeterHit pSimHit = pMCRecoCaloAssociationCollection.at(ic).getSim();
            for (int iCont = 0, iEnd = pSimHit.contributions_size(); iCont < iEnd; ++iCont) {
              edm4hep::ConstCaloHitContribution conb = pSimHit.getContributions(iCont);
              const edm4hep::ConstMCParticle ipa = conb.getParticle();
              float ien = conb.getEnergy();

              if (m_id_pMC_map->find(ipa.id()) == m_id_pMC_map->end())
                continue;

              const edm4hep::MCParticle* pMCParticle = (*m_id_pMC_map)[ipa.id()];
              mcParticleToEnergyWeightMap[pMCParticle] += ien;
            }
          }

          for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter = mcParticleToEnergyWeightMap.begin(),
                                                           mcParticleIterEnd = mcParticleToEnergyWeightMap.end();
               mcParticleIter != mcParticleIterEnd; ++mcParticleIter) {
            PANDORA_THROW_RESULT_IF(
                pandora::STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetCaloHitToMCParticleRelationship(m_pandora, &calorimeterHitVector.at(i_calo),
                                                               mcParticleIter->first, mcParticleIter->second));
          }
        } catch (const pandora::StatusCodeException& statusCodeException) {
          streamlog_out(ERROR) << "Failed to extract calo hit to mc particle relationship: "
                               << statusCodeException.ToString() << std::endl;
        } catch (const std::exception& exception) {
          streamlog_out(WARNING) << "Failed to extract calo hit to mc particle relationship: " << exception.what()
                                 << std::endl;
        }
      }
    } catch (const std::exception& exception) {
      streamlog_out(DEBUG5) << "Failed to extract calo hit to mc particle relationships collection: " << *iter << ", "
                            << exception.what() << std::endl;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::Settings::Settings() {}

