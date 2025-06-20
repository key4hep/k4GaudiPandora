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
#include "DDMCParticleCreator.h"

#include <edm4hep/CaloHitContribution.h>
#include <edm4hep/CaloHitMCParticleLinkCollection.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/Track.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>

#include <Objects/Helix.h>

#include <cmath>
#include <limits>
#include <stdexcept>

// forward declarations. See in DDPandoraPFANewProcessor.cc
double getFieldFromCompact();

DDMCParticleCreator::DDMCParticleCreator(const Settings& settings, pandora::Pandora& pandora)
    : m_settings(settings), m_pandora(pandora), m_bField(getFieldFromCompact()) {}

DDMCParticleCreator::~DDMCParticleCreator() {}

pandora::StatusCode
DDMCParticleCreator::CreateMCParticles(const std::vector<edm4hep::MCParticle>& mcParticleCollections) const {

  for (const auto& mcParticle : mcParticleCollections) {
    PandoraApi::MCParticle::Parameters mcParticleParameters;
    mcParticleParameters.m_energy = mcParticle.getEnergy();
    mcParticleParameters.m_particleId = mcParticle.getPDG();
    mcParticleParameters.m_mcParticleType = pandora::MC_3D;
    mcParticleParameters.m_pParentAddress = &mcParticle;
    mcParticleParameters.m_momentum =
        pandora::CartesianVector(mcParticle.getMomentum().x, mcParticle.getMomentum().y, mcParticle.getMomentum().z);
    mcParticleParameters.m_vertex =
        pandora::CartesianVector(mcParticle.getVertex().x, mcParticle.getVertex().y, mcParticle.getVertex().z);
    mcParticleParameters.m_endpoint =
        pandora::CartesianVector(mcParticle.getEndpoint().x, mcParticle.getEndpoint().y, mcParticle.getEndpoint().z);

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters))

    // Create parent-daughter relationships
    for (const auto& daughter : mcParticle.getDaughters()) {
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::SetMCParentDaughterRelationship(m_pandora, &mcParticle, &daughter))
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

pandora::StatusCode
DDMCParticleCreator::CreateTrackToMCParticleRelationships(const MCPCollectionVector& mcParticleCollections,
                                                          const TrackMCLinkCollectionVector& trackRelCollections,
                                                          const TrackVector& trackVector) const {
  for (const auto& trackRelCollection : trackRelCollections) {
    for (const auto& pTrack : trackVector) {
      try {
        // Get reconstructed momentum at dca
        const auto& trackState = pTrack.getTrackStates(0);
        const pandora::Helix helixFit(trackState.phi, trackState.D0, trackState.Z0, trackState.omega,
                                      trackState.tanLambda, m_bField);
        const float recoMomentum = helixFit.GetMomentum().GetMagnitude();

        // APS: I am not sure why we chose the best MCParticle per trackRelationCollection, but this was like this from
        // the start, and we only have one track collection, so it doesn't matter...

        // Use momentum magnitude to identify best MC particle
        const edm4hep::MCParticle* pBestMCParticle = nullptr;
        float bestDeltaMomentum = std::numeric_limits<float>::max();

        for (const auto&& trackMCRel : *trackRelCollection) {
          // if the track relation does not match the track, continue
          if (trackMCRel.get<edm4hep::Track>().id() != pTrack.id())
            continue;

          auto const& mcParticle = trackMCRel.get<edm4hep::MCParticle>();

          const float trueMomentum = pandora::CartesianVector(mcParticle.getMomentum()[0], mcParticle.getMomentum()[1],
                                                              mcParticle.getMomentum()[2])
                                         .GetMagnitude();
          const float deltaMomentum = std::fabs(recoMomentum - trueMomentum);

          if (deltaMomentum < bestDeltaMomentum) {
            pBestMCParticle = &mcParticle;
            bestDeltaMomentum = deltaMomentum;
          }
        }

        if (pBestMCParticle == nullptr) {
          // streamlog_out(WARNING) << "No suitable MC particle found for track association." << std::endl;
          continue;
        }

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                PandoraApi::SetTrackToMCParticleRelationship(m_pandora, &pTrack, pBestMCParticle))

      } catch (const pandora::StatusCodeException& statusCodeException) {
        // m_algorithm.error() << "Failed to extract track to MC particle relationship: "
        //                      << statusCodeException.ToString() << std::endl;
      } catch (const std::exception& exception) {
        // streamlog_out(WARNING) << "Exception encountered: " << exception.what() << std::endl;
      }
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateCaloHitToMCParticleRelationships(
    const CaloHitSimCaloHitLinkCollectionVector& caloRelCollections, const HitVector& calorimeterHitVector) const {
  using MCParticleToEnergyWeightMap = std::map<const edm4hep::MCParticle*, float>;
  MCParticleToEnergyWeightMap mcParticleToEnergyWeightMap;

  for (const auto& hitRelCollection : caloRelCollections) {
    for (const auto& caloHit : calorimeterHitVector) {
      try {
        for (const auto&& caloHitLink : *hitRelCollection) {
          try {
            mcParticleToEnergyWeightMap.clear();

            if (caloHitLink.get<edm4hep::CalorimeterHit>().id() != caloHit.id())
              continue;

            const auto& simHit = caloHitLink.get<edm4hep::SimCalorimeterHit>();
            for (const auto& conb : simHit.getContributions()) {
              const auto& ipa = conb.getParticle();
              float ien = conb.getEnergy();

              mcParticleToEnergyWeightMap[&ipa] += ien;
            }

            for (const auto& mcPToE : mcParticleToEnergyWeightMap) {
              PANDORA_THROW_RESULT_IF(
                  pandora::STATUS_CODE_SUCCESS, !=,
                  PandoraApi::SetCaloHitToMCParticleRelationship(m_pandora, &caloHit, mcPToE.first, mcPToE.second));
            }
          } catch (const pandora::StatusCodeException& statusCodeException) {
            // m_algorithm.error() << "Failed to extract calo hit to mc particle relationship: "
            //                      << statusCodeException.ToString() << std::endl;
          } catch (const std::exception& exception) {
            // streamlog_out(WARNING) << "Failed to extract calo hit to mc particle relationship: " << exception.what()
            //                        << std::endl;
          }
        }
      } catch (const std::exception& exception) {
        // streamlog_out(DEBUG5) << "Failed to extract calo hit to mc particle relationships collection: " << *iter <<
        // ", "
        //                       << exception.what() << std::endl;
      }
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::Settings::Settings()
    : m_mcParticleCollections(StringVector()), m_caloHitRelationCollections(StringVector()),
      m_trackRelationCollections(StringVector()) {}
