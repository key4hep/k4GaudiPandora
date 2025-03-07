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
 *  @file   DDMarlinPandora/include/DDPandoraPFANewAlgorithm.h
 *
 *  @brief  Header file for the pandora pfa new algorithm class.
 *
 *  $Log: $
 */

#ifndef DDPANDORAPFANEWALGORITHM_H
#define DDPANDORAPFANEWALGORITHM_H 1

#include "DDCaloHitCreator.h"
#include "DDGeometryCreator.h"
#include "DDMCParticleCreator.h"
#include "DDPfoCreator.h"
#include "DDTrackCreatorBase.h"

// k4FWCore
#include <k4FWCore/DataHandle.h>
#include <k4FWCore/BaseClass.h>
#include <k4FWCore/Transformer.h>
#include <Gaudi/Property.h>

#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"

namespace pandora {
  class Pandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDPandoraPFANewAlgorithm class
 */
struct DDPandoraPFANewAlgorithm final : 
  k4FWCore::MultiTransformer<std::tuple<edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection, edm4hep::VertexCollection>(
	const std::vector<const edm4hep::MCParticle*>&,
   const std::vector<const edm4hep::VertexCollection*>&,
   const std::vector<const edm4hep::VertexCollection*>&,
   const std::vector<const edm4hep::VertexCollection*>&,
   const std::vector<const edm4hep::TrackCollection*>&,
   const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,
   const std::vector<const edm4hep::CalorimeterHitCollection*>&,
   const std::vector<const edm4hep::CalorimeterHitCollection*>&,
   const std::vector<const edm4hep::CalorimeterHitCollection*>&,
   const std::vector<const edm4hep::CalorimeterHitCollection*>&,
   const std::vector<const edm4hep::CalorimeterHitCollection*>&,
   const std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*>&)> {
public:
  typedef std::vector<float>       FloatVector;
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

    std::string m_pandoraSettingsXmlFile = "";  ///< The pandora settings xml file

    float m_innerBField      = 0.0;    ///< The bfield in the main tracker, ecal and hcal, units Tesla
    float m_muonBarrelBField = 0.0;    ///< The bfield in the muon barrel, units Tesla
    float m_muonEndCapBField = 0.0;    ///< The bfield in the muon endcap, units Tesla
    bool  m_useDD4hepField   = false;  ///< Whether to use the DD4hep field map instead of the values above

    FloatVector m_inputEnergyCorrectionPoints{};   ///< The input energy points for non-linearity energy correction
    FloatVector m_outputEnergyCorrectionPoints{};  ///< The output energy points for non-linearity energy correction

    // Software compensation parameters
    FloatVector m_softCompParameters{};
    FloatVector m_softCompEnergyDensityBins{};
    float       m_energyDensityFinalBin           = 0.0;
    float       m_maxClusterEnergyToApplySoftComp = 100.0;
    float       m_minCleanHitEnergy               = 0.5;
    float       m_minCleanHitEnergyFraction       = 0.01;
    float       m_minCleanCorrectedHitEnergy      = 0.1;

    ///ADDED BY NIKIFOROS
    //Detector names not needed anymore, accessed by det type flags
    std::string m_trackCreatorName = "";  ///< The name of the DDTrackCreator implementation to use
  };

  /**
     *  @brief  Default constructor
     */
  DDPandoraPFANewAlgorithm(const std::string& name, ISvcLocator* svcLoc);

  /**
     *  @brief  Initialize, called at startup
     */
  StatusCode initialize();


  /**
   *  @brief operator, the workhorse of the algorithm
   * 
   *  @param  MCParticleCollections Collection of MCParticles
   *  @param  kinkCollections the Vertex collections of kinks
   *  @param  prongCollections the Vertex collections of prongs and splits
   *  @param  v0Collections the Vertex collections of V0s
   *  @param  trackerHitLinkCollections the associations between trackerHits and simTrackerHits
   *  @param  trackCollections collections of tracks
   *  @param  eCalCollections  CalorimeterHit Collection for the ECal
   *  @param  hCalCollections  CalorimeterHit Collection for the HCal
   *  @param  mCalCollections  CalorimeterHit Collection for the Muon Calo
   *  @param  lCalCollections  CalorimeterHit Collection for the LCal
   *  @param  lhCalCollections CalorimeterHit Collection for the HLCal
   *  @param  caloLinkCollections the associations between CalorimeterHits and MCParticles
   * 
   *  @return tuple of reconstructed: (Clusters, RecoParticles, Verticies)
   *  
   */
  std::tuple<edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection, edm4hep::VertexCollection> operator()(
   const std::vector<const edm4hep::MCParticle*>& MCParticleCollections,
   const std::vector<const edm4hep::VertexCollection*>& kinkCollections,
   const std::vector<const edm4hep::VertexCollection*>& prongCollections,
   const std::vector<const edm4hep::VertexCollection*>& v0Collections,
   const std::vector<const edm4hep::TrackCollection*>& trackCollections,
   const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& trackerHitLinkCollections,
   const std::vector<const edm4hep::CalorimeterHitCollection*>& eCalCollections,
   const std::vector<const edm4hep::CalorimeterHitCollection*>& hCalCollections,
   const std::vector<const edm4hep::CalorimeterHitCollection*>& mCalCollections,
   const std::vector<const edm4hep::CalorimeterHitCollection*>& lCalCollections,
   const std::vector<const edm4hep::CalorimeterHitCollection*>& lhCalCollections,
   const std::vector<const edm4hep::CaloHitMCParticleLinkCollection*>& caloLinkCollections
 ) const override;



  /**
     *  @brief  End, called at shutdown
     */
  StatusCode finalize();

  /**
     *  @brief  Get address of the pandora instance
     *
     *  @return address of the pandora instance
     */
  const pandora::Pandora* GetPandora() const;

  /**
     *  @brief  Get address of the current lcio event
     *
     *  @param  pPandora address of the relevant pandora instance
     *
     *  @return address of the current lcio event
     */
  //static const EVENT::LCEvent* GetCurrentEvent(const pandora::Pandora* const pPandora);

private:
  /**
     *  @brief  Register user algorithm factories, energy correction functions and particle id functions,
     *          insert user code here
     */
  pandora::StatusCode RegisterUserComponents() const;

  /**
     *  @brief  Process steering file parameters, insert user code here
     */
  void ProcessSteeringFile();

  /**
     *  @brief  Copy some steering parameters between settings objects
     */
  void FinaliseSteeringParameters();

  /**
     *  @brief  Reset the pandora pfa new processor
     */
  void Reset();

  pandora::Pandora*    m_pPandora             = NULL;  ///< Address of the pandora instance
  DDCaloHitCreator*    m_pCaloHitCreator      = NULL;  ///< The calo hit creator
  DDGeometryCreator*   m_pGeometryCreator     = NULL;  ///< The geometry creator
  DDTrackCreatorBase*  m_pTrackCreator        = NULL;  ///< The track creator
  DDMCParticleCreator* m_pDDMCParticleCreator = NULL;  ///< The mc particle creator
  DDPfoCreator*        m_pDDPfoCreator        = NULL;  ///< The pfo creator

  Settings                      m_settings{};                   ///< The settings for the pandora pfa new processor
  DDCaloHitCreator::Settings    m_caloHitCreatorSettings{};     ///< The calo hit creator settings
  DDGeometryCreator::Settings   m_geometryCreatorSettings{};    ///< The geometry creator settings
  DDMCParticleCreator::Settings m_mcParticleCreatorSettings{};  ///< The mc particle creator settings
  DDTrackCreatorBase::Settings  m_trackCreatorSettings{};       ///< The track creator settings
  DDPfoCreator::Settings        m_pfoCreatorSettings{};         ///< The pfo creator settings

  //typedef std::map<const pandora::Pandora*, EVENT::LCEvent*> PandoraToLCEventMap;
  //static PandoraToLCEventMap                                 m_pandoraToLCEventMap;  ///< The pandora to lc event map
};

//------------------------------------------------------------------------------------------------------------------------------------------

#endif  // #ifndef DDPANDORAPFANEWALGORITHM_H
