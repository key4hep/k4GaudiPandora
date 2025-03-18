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

#ifndef DDMCPARTICLECREATOR_H
#define DDMCPARTICLECREATOR_H 1

#include "Api/PandoraApi.h"

#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>

#include "DDCaloHitCreator.h"
#include "DDTrackCreatorBase.h"
/**
 *  @brief  DDMCParticleCreator class
 */
class CollectionMaps;
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

    float        m_bField;                        ///< m_bField
  };

  /**
     *  @brief  Constructor
     *
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
  DDMCParticleCreator(const Settings& settings, const pandora::Pandora* const pPandora, IMessageSvc* msgSvc);

  /**
     *  @brief  Destructor
     */
  ~DDMCParticleCreator();

  /**
     *  @brief  Create MCParticles
     *
     *  @param  pLCEvent the lcio event
     */
  pandora::StatusCode CreateMCParticles(const std::vector<const edm4hep::MCParticleCollection*>& MCParticleCollections);

  /**
     *  @brief  Create Track to mc particle relationships
     *
     */
  pandora::StatusCode CreateTrackToMCParticleRelationships(const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& linkCollections,
                                                           const TrackVector&    trackVector,
                                                           const std::vector<const edm4hep::TrackCollection*>& trackCollections) const;

  /**
     *  @brief  Create calo hit to mc particle relationships
     *
     */
  pandora::StatusCode CreateCaloHitToMCParticleRelationships(const std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*>& linkCollections,
                                                             const CalorimeterHitVector& calorimeterHitVector,
                                                             const std::vector<const edm4hep::CalorimeterHitCollection*>& eCalCollections,
                                                             const std::vector<const edm4hep::CalorimeterHitCollection*>& hCalCollections,
                                                             const std::vector<const edm4hep::CalorimeterHitCollection*>& mCalCollections,
                                                             const std::vector<const edm4hep::CalorimeterHitCollection*>& lCalCollections,
                                                             const std::vector<const edm4hep::CalorimeterHitCollection*>& lhCalCollections) const;
  /**
     *  @brief  Reset the MCP creator
     */
  void Reset();

private:
  const Settings          m_settings;  ///< The mc particle creator settings
  const pandora::Pandora& m_pandora;   ///< Reference to the pandora object to create the mc particles
  const float             m_bField;    ///< The bfield  

  IMessageSvc* m_msgSvc;
};

#endif  // #ifndef DDMCPARTICLECREATOR_H
