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
 *  @file   DDMarlinPandora/include/DDMCParticleCreator.h
 *
 *  @brief  Header file for the mc particle creator class.
 *
 *  $Log: $
 */

/**
 *  @brief  DDMCParticleCreator class
 */
#ifndef DDMCPARTICLECREATOR_H
#define DDMCPARTICLECREATOR_H

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "DDCaloHitCreator.h"
#include "DDTrackCreatorBase.h"

#include "Api/PandoraApi.h"
#include <vector>
#include <map>
#include <string>

class CollectionMaps;

/**
 *  @brief  DDMCParticleCreator class
 */
class DDMCParticleCreator {
public:
  typedef std::vector<std::string> StringVector;

  /**
   *  @brief  Settings class
   */
  class Settings {
  public:
    /**
     *  @brief  Default constructor
     */
    Settings();

    StringVector m_mcParticleCollections;         ///< The mc particle collections
    StringVector m_caloHitRelationCollections;    ///< The SimCaloHit to CaloHit particle relations
    StringVector m_trackRelationCollections;      ///< The SimTrackerHit to TrackerHit particle relations
    float        m_bField;                        ///< Magnetic field strength
  };

  /**
   *  @brief  Constructor
   *
   *  @param  settings the creator settings
   *  @param  pPandora address of the relevant pandora instance
   */
  DDMCParticleCreator(const Settings& settings, const pandora::Pandora* const pPandora);

  /**
   *  @brief  Destructor
   */
  ~DDMCParticleCreator();

  /**
   *  @brief  Create MCParticles
   *
   *  @param  collectionMaps The collection map containing input data
   */
  pandora::StatusCode CreateMCParticles(const CollectionMaps& collectionMaps) const;

  /**
   *  @brief  Create Track to MCParticle relationships
   *
   *  @param  collectionMaps The collection map containing input data
   *  @param  trackVector Vector of reconstructed tracks
   */
  pandora::StatusCode CreateTrackToMCParticleRelationships(const CollectionMaps& collectionMaps,
                                                           const std::vector<edm4hep::Track>& trackVector) const;

  /**
   *  @brief  Create CaloHit to MCParticle relationships
   *
   *  @param  collectionMaps The collection map containing input data
   *  @param  calorimeterHitVector Vector of calorimeter hits
   */
  pandora::StatusCode CreateCaloHitToMCParticleRelationships(const CollectionMaps& collectionMaps,
                                                             const std::vector<edm4hep::CalorimeterHit>& calorimeterHitVector) const;

private:
  const Settings          m_settings;  ///< The mc particle creator settings
  const pandora::Pandora& m_pandora;   ///< Reference to the pandora object to create the mc particles
  const float             m_bField;    ///< The magnetic field strength
  std::map<unsigned int, const edm4hep::MCParticle*>* m_id_pMC_map;
};

#endif  // #ifndef DDMCPARTICLECREATOR_H

