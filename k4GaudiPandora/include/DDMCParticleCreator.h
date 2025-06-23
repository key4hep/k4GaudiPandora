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

#ifndef K4GAUDIPANDORA_DDMCPARTICLECREATOR_H
#define K4GAUDIPANDORA_DDMCPARTICLECREATOR_H

#include <Api/PandoraApi.h>

#include <edm4hep/CaloHitMCParticleLinkCollection.h>
#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>

#include "Gaudi/Algorithm.h"

#include <string>
#include <vector>

typedef std::vector<const edm4hep::MCParticleCollection*> MCPCollectionVector;
typedef std::vector<const edm4hep::TrackMCParticleLinkCollection*> TrackMCLinkCollectionVector;
typedef std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*> CaloHitSimCaloHitLinkCollectionVector;
typedef std::vector<edm4hep::Track> TrackVector;
typedef std::vector<edm4hep::CalorimeterHit> HitVector;

class DDMCParticleCreator {
public:
  typedef std::vector<std::string> StringVector;

  class Settings {
  public:
    Settings();

    StringVector m_mcParticleCollections;      ///< The mc particle collections
    StringVector m_caloHitRelationCollections; ///< The SimCaloHit to CaloHit particle relations
    StringVector m_trackRelationCollections;   ///< The Track to MCParticle relation collections
    float m_bField;                            ///< Magnetic field strength
  };

  /**
   *  @brief  Constructor
   *
   *  @param  settings the creator settings
   *  @param  pandora reference to the relevant pandora instance
   *  @param  algorithm pointer to the Gaudi algorithm instance
   */
  DDMCParticleCreator(const Settings& settings, pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm);

  pandora::StatusCode CreateMCParticles(const std::vector<edm4hep::MCParticle>& mcParticleCollections) const;

  /**
   *  @brief  Create Track to MCParticle relationships
   *
   *  @param  collectionMaps The collection map containing input data
   *  @param  trackVector Vector of reconstructed tracks
   */
  pandora::StatusCode CreateTrackToMCParticleRelationships(const TrackMCLinkCollectionVector& trackRelCollections,
                                                           const TrackVector& trackVector) const;

  /**
   *  @brief  Create CaloHit to MCParticle relationships
   *
   *  @param  collectionMaps The collection map containing input data
   *  @param  calorimeterHitVector Vector of calorimeter hits
   */
  pandora::StatusCode
  CreateCaloHitToMCParticleRelationships(const CaloHitSimCaloHitLinkCollectionVector& caloRelCollections,
                                         const HitVector& calorimeterHitVector) const;

private:
  const Settings m_settings;           ///< The mc particle creator settings
  pandora::Pandora& m_pandora;         ///< Reference to the pandora object to create the mc particles
  const float m_bField;                ///< The magnetic field strength
  const Gaudi::Algorithm& m_algorithm; ///< Reference to the Gaudi algorithm for logging
};

#endif // #ifndef K4GAUDIPANDORA_DDMCPARTICLECREATOR_H
