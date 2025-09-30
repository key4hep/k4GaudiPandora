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

#ifndef DDPANDORAPFANEWALGORITHM_H
#define DDPANDORAPFANEWALGORITHM_H 1

#include "DDCaloHitCreator.h"
#include "DDGeometryCreator.h"
#include "DDMCParticleCreator.h"
#include "DDPfoCreator.h"
#include "DDTrackCreatorBase.h"

#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/ClusterCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>
#include <edm4hep/VertexCollection.h>

#include <Gaudi/Property.h>
#include <k4FWCore/Transformer.h>
#include <k4Interface/IGeoSvc.h>

#include <LCPlugins/LCSoftwareCompensation.h>

#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace pandora {
class Pandora;
}

dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag = 0);

struct DDPandoraPFANewAlgorithm final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection,
                                            edm4hep::VertexCollection>(
          const std::vector<const edm4hep::MCParticleCollection*>&,
          const std::vector<const edm4hep::VertexCollection*>&, const std::vector<const edm4hep::VertexCollection*>&,
          const std::vector<const edm4hep::VertexCollection*>&, const std::vector<const edm4hep::VertexCollection*>&,
          const std::vector<const edm4hep::TrackCollection*>&,
          const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,
          const std::vector<const edm4hep::CalorimeterHitCollection*>&,
          const std::vector<const edm4hep::CalorimeterHitCollection*>&,
          const std::vector<const edm4hep::CalorimeterHitCollection*>&,
          const std::vector<const edm4hep::CalorimeterHitCollection*>&,
          const std::vector<const edm4hep::CalorimeterHitCollection*>&,
          const std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*>&)> {
public:
  class Settings {
  public:
    Settings();

    std::string m_pandoraSettingsXmlFile = ""; ///< The pandora settings xml file

    float m_innerBField = 0.0;      ///< The bfield in the main tracker, ecal and hcal, units Tesla
    float m_muonBarrelBField = 0.0; ///< The bfield in the muon barrel, units Tesla
    float m_muonEndCapBField = 0.0; ///< The bfield in the muon endcap, units Tesla
    bool m_useDD4hepField = false;  ///< Whether to use the DD4hep field map instead of the values above

    std::vector<float> m_inputEnergyCorrectionPoints{}; ///< The input energy points for non-linearity energy correction
    std::vector<float>
        m_outputEnergyCorrectionPoints{}; ///< The output energy points for non-linearity energy correction
    std::vector<float> m_ecalInputEnergyCorrectionPoints{};  ///< The input energy points for non-linearity energy
                                                             ///< correction in the ECAL
    std::vector<float> m_ecalOutputEnergyCorrectionPoints{}; ///< The input energy points for non-linearity energy
                                                             ///< correction in the ECAL

    // Software compensation parameters
    std::vector<float> m_softCompParameters{};
    std::vector<float> m_softCompEnergyDensityBins{};
    float m_energyDensityFinalBin = 0.0;
    float m_maxClusterEnergyToApplySoftComp = 100.0;
    float m_minCleanHitEnergy = 0.5;
    float m_minCleanHitEnergyFraction = 0.01;
    float m_minCleanCorrectedHitEnergy = 0.1;

    /// ADDED BY NIKIFOROS
    // Detector names not needed anymore, accessed by det type flags
    std::string m_trackCreatorName = ""; ///< The name of the DDTrackCreator implementation to use

    // Detector name (used by ALLEGRO)
    std::string m_detectorName = "";
  };

  DDPandoraPFANewAlgorithm(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  /**
   *  @brief operator, the workhorse of the algorithm
   *
   *  @param  MCParticleCollections Collection of MCParticles
   *  @param  kinkCollections the Vertex collections of kinks
   *  @param  prongCollections the Vertex collections of prongs
   *  @param  splitCollections the Vertex collections of splits
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
  std::tuple<edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection, edm4hep::VertexCollection>
  operator()(const std::vector<const edm4hep::MCParticleCollection*>& MCParticleCollections,
             const std::vector<const edm4hep::VertexCollection*>& kinkCollections,
             const std::vector<const edm4hep::VertexCollection*>& prongCollections,
             const std::vector<const edm4hep::VertexCollection*>& splitCollections,
             const std::vector<const edm4hep::VertexCollection*>& v0Collections,
             const std::vector<const edm4hep::TrackCollection*>& trackCollections,
             const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& trackerHitLinkCollections,
             const std::vector<const edm4hep::CalorimeterHitCollection*>& eCalCollections,
             const std::vector<const edm4hep::CalorimeterHitCollection*>& hCalCollections,
             const std::vector<const edm4hep::CalorimeterHitCollection*>& mCalCollections,
             const std::vector<const edm4hep::CalorimeterHitCollection*>& lCalCollections,
             const std::vector<const edm4hep::CalorimeterHitCollection*>& lhCalCollections,
             const std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*>& caloLinkCollections) const override;

  const pandora::Pandora* GetPandora() const;

private:
  /**
   *  @brief  Register user algorithm factories, energy correction functions and particle id functions,
   *          insert user code here
   */
  pandora::StatusCode registerUserComponents() const;

  /**
   *  @brief  Copy some steering parameters between settings objects
   */
  void finaliseSteeringParameters();

  /**
   *  @brief  Reset the pandora pfa new processor
   */
  void reset() const;

  pandora::Pandora m_pPandora;                                 ///< Pandora instance
  std::unique_ptr<DDCaloHitCreator> m_caloHitCreator;          ///< The calo hit creator
  std::unique_ptr<DDGeometryCreator> m_geometryCreator;        ///< The geometry creator
  std::unique_ptr<DDTrackCreatorBase> m_pTrackCreator;         ///< The track creator
  std::unique_ptr<DDMCParticleCreator> m_pDDMCParticleCreator; ///< The mc particle creator
  std::unique_ptr<DDPfoCreator> m_pfoCreator;                  ///< The pfo creator
  SmartIF<IGeoSvc> m_geoSvc;                                   ///< The GeoSvc

  Settings m_settings{};                                       ///< The settings for the pandora pfa new processor
  DDCaloHitCreator::Settings m_caloHitCreatorSettings{};       ///< The calo hit creator settings
  DDGeometryCreator::Settings m_geometryCreatorSettings{};     ///< The geometry creator settings
  DDMCParticleCreator::Settings m_mcParticleCreatorSettings{}; ///< The mc particle creator settings
  DDTrackCreatorBase::Settings m_trackCreatorSettings{};       ///< The track creator settings
  DDPfoCreator::Settings m_pfoCreatorSettings{};               ///< The pfo creator settings

  Gaudi::Property<std::string> m_pandoraSettingsXmlFile{this, "PandoraSettingsXmlFile", "",
                                                        "The pandora settings xml file"};
  Gaudi::Property<bool> m_createGaps{this, "CreateGaps", true,
                                     "Decides whether to create gaps in the geometry (ILD-specific)"};
  Gaudi::Property<std::string> m_startVertexAlgName{this, "StartVertexAlgorithmName", "PandoraPFANew",
                                                    "The algorithm name for filling start vertex"};
  // Energy resolution parameters
  Gaudi::Property<float> m_emStochasticTerm{this, "EMStochasticTerm", 0.17f, "The stochastic term for EM shower"};
  Gaudi::Property<float> m_hadStochasticTerm{this, "HadStochasticTerm", 0.6f,
                                             "The stochastic term for Hadronic shower"};
  Gaudi::Property<float> m_emConstantTerm{this, "EMConstantTerm", 0.01f, "The constant term for EM shower"};
  Gaudi::Property<float> m_hadConstantTerm{this, "HadConstantTerm", 0.03f, "The constant term for Hadronic shower"};
  // Calibration constants
  Gaudi::Property<float> m_eCalToMip{this, "ECalToMipCalibration", 1.0f,
                                     "The calibration from deposited ECal energy to mip"};
  Gaudi::Property<float> m_hCalToMip{this, "HCalToMipCalibration", 1.0f,
                                     "The calibration from deposited HCal energy to mip"};
  Gaudi::Property<float> m_eCalMipThreshold{this, "ECalMipThreshold", 0.0f,
                                            "Threshold for creating calo hits in the ECal, units mip"};
  Gaudi::Property<float> m_muonToMip{this, "MuonToMipCalibration", 1.0f,
                                     "The calibration from deposited Muon energy to mip"};
  Gaudi::Property<float> m_hCalMipThreshold{this, "HCalMipThreshold", 0.0f,
                                            "Threshold for creating calo hits in the HCal, units mip"};
  Gaudi::Property<float> m_eCalToEMGeV{this, "ECalToEMGeVCalibration", 1.0f,
                                       "The calibration from deposited ECal energy to EM energy"};
  Gaudi::Property<float> m_hCalToEMGeV{this, "HCalToEMGeVCalibration", 1.0f,
                                       "The calibration from deposited HCal energy to EM energy"};
  Gaudi::Property<float> m_eCalToHadGeVEndCap{this, "ECalToHadGeVCalibrationEndCap", 1.0f,
                                              "The calibration from deposited ECal energy to hadronic energy"};
  Gaudi::Property<float> m_eCalToHadGeVBarrel{this, "ECalToHadGeVCalibrationBarrel", 1.0f,
                                              "The calibration from deposited ECal energy to hadronic energy"};
  Gaudi::Property<float> m_hCalToHadGeV{this, "HCalToHadGeVCalibration", 1.0f,
                                        "The calibration from deposited HCal energy to hadronic energy"};
  Gaudi::Property<int> m_muonDigitalHits{this, "DigitalMuonHits", 1, "Treat muon hits as digital"};
  Gaudi::Property<float> m_muonHitEnergy{this, "MuonHitEnergy", 0.5f,
                                         "The energy for a digital muon calorimeter hit, units GeV"};
  Gaudi::Property<float> m_maxHCalHitHadronicEnergy{this, "MaxHCalHitHadronicEnergy", 10000.0f,
                                                    "The maximum hadronic energy allowed for a single hcal hit"};
  Gaudi::Property<int> m_nOuterSamplingLayers{this, "NOuterSamplingLayers", 3,
                                              "Number of layers from edge for hit to be flagged as an outer layer hit"};
  Gaudi::Property<float> m_layersFromEdgeMaxRearDistance{
      this, "LayersFromEdgeMaxRearDistance", 250.0f,
      "Maximum number of layers from candidate outer layer hit to rear of detector"};
  // B-field parameters
  Gaudi::Property<float> m_muonBarrelBField{this, "MuonBarrelBField", -1.5f,
                                            "The bfield in the muon barrel, units Tesla"};
  Gaudi::Property<float> m_muonEndCapBField{this, "MuonEndCapBField", 0.01f,
                                            "The bfield in the muon endcap, units Tesla"};
  Gaudi::Property<bool> m_useDD4hepField{this, "UseDD4hepField", false, "Whether to use the BField map from DD4hep"};
  // Track relationship parameters
  Gaudi::Property<int> m_shouldFormTrackRelationships{
      this, "ShouldFormTrackRelationships", 1, "Whether to form pandora track relationships using v0 and kink info"};
  // Initial track hit specifications
  Gaudi::Property<int> m_minTrackHits{this, "MinTrackHits", 5, "Track quality cut: the minimum number of track hits"};
  Gaudi::Property<int> m_minFtdTrackHits{this, "MinFtdTrackHits", 0,
                                         "Track quality cut: the minimum number of ftd track hits for ftd only tracks"};
  Gaudi::Property<int> m_maxTrackHits{this, "MaxTrackHits", 5000,
                                      "Track quality cut: the maximum number of track hits"};
  // Track PFO usage parameters
  Gaudi::Property<float> m_d0TrackCut{this, "D0TrackCut", 50.0f,
                                      "Track d0 cut used to determine whether track can be used to form pfo"};
  Gaudi::Property<float> m_z0TrackCut{this, "Z0TrackCut", 50.0f,
                                      "Track z0 cut used to determine whether track can be used to form pfo"};
  Gaudi::Property<int> m_usingNonVertexTracks{this, "UseNonVertexTracks", 1,
                                              "Whether can form pfos from tracks that don't start at vertex"};
  Gaudi::Property<int> m_usingUnmatchedNonVertexTracks{
      this, "UseUnmatchedNonVertexTracks", 0, "Whether can form pfos from unmatched tracks that don't start at vertex"};
  Gaudi::Property<int> m_usingUnmatchedVertexTracks{this, "UseUnmatchedVertexTracks", 1,
                                                    "Whether can form pfos from unmatched tracks that start at vertex"};
  Gaudi::Property<float> m_unmatchedVertexTrackMaxEnergy{this, "UnmatchedVertexTrackMaxEnergy", 5.0f,
                                                         "Maximum energy for unmatched vertex track"};
  Gaudi::Property<float> m_d0UnmatchedVertexTrackCut{
      this, "D0UnmatchedVertexTrackCut", 5.0f, "d0 cut used to determine whether unmatched vertex track can form pfo"};
  Gaudi::Property<float> m_z0UnmatchedVertexTrackCut{
      this, "Z0UnmatchedVertexTrackCut", 5.0f, "z0 cut used to determine whether unmatched vertex track can form pfo"};
  Gaudi::Property<float> m_zCutForNonVertexTracks{
      this, "ZCutForNonVertexTracks", 250.0f, "Non vtx track z cut to determine whether track can be used to form pfo"};
  // Track "reaches ecal" parameters
  Gaudi::Property<int> m_reachesECalNBarrelTrackerHits{
      this, "ReachesECalNBarrelTrackerHits", 11,
      "Minimum number of BarrelTracker hits to consider track as reaching ecal"};
  Gaudi::Property<int> m_reachesECalNFtdHits{this, "ReachesECalNFtdHits", 4,
                                             "Minimum number of FTD hits to consider track as reaching ecal"};
  Gaudi::Property<float> m_reachesECalBarrelTrackerOuterDistance{
      this, "ReachesECalBarrelTrackerOuterDistance", -100.0f,
      "Max distance from track to BarrelTracker r max to determine whether track reaches ecal"};
  Gaudi::Property<int> m_reachesECalMinFtdLayer{this, "ReachesECalMinFtdLayer", 9,
                                                "Min FTD layer for track to be considered to have reached ecal"};
  Gaudi::Property<float> m_reachesECalBarrelTrackerZMaxDistance{
      this, "ReachesECalBarrelTrackerZMaxDistance", -50.0f,
      "Max distance from track to BarrelTracker z max to determine whether track reaches ecal"};
  Gaudi::Property<float> m_reachesECalFtdZMaxDistance{
      this, "ReachesECalFtdZMaxDistance", 1.0f, "Max distance from track hit to FTD z position to identify FTD hits"};
  Gaudi::Property<float> m_curvatureToMomentumFactor{this, "CurvatureToMomentumFactor", 0.3f / 2000.0f,
                                                     "Constant relating track curvature in B field to momentum"};
  Gaudi::Property<float> m_minTrackECalDistanceFromIp{
      this, "MinTrackECalDistanceFromIp", 100.0f,
      "Sanity check on separation between IP and track projected ECal position"};
  // Final track quality parameters
  Gaudi::Property<float> m_maxTrackSigmaPOverP{this, "MaxTrackSigmaPOverP", 0.15f,
                                               "Cut on fractional track momentum error"};
  Gaudi::Property<float> m_minMomentumForTrackHitChecks{
      this, "MinMomentumForTrackHitChecks", 1.0f,
      "Min track momentum required to perform final quality checks on number of hits"};
  Gaudi::Property<float> m_minBarrelTrackerHitFractionOfExpected{
      this, "MinBarrelTrackerHitFractionOfExpected", 0.20f, "Cut on fraction of expected number of BarrelTracker hits"};
  Gaudi::Property<int> m_minFtdHitsForBarrelTrackerHitFraction{
      this, "MinFtdHitsForBarrelTrackerHitFraction", 2,
      "Cut on minimum number of FTD hits for BarrelTracker hit fraction to be applied"};
  Gaudi::Property<float> m_maxBarrelTrackerInnerRDistance{
      this, "MaxBarrelTrackerInnerRDistance", 50.0f,
      "Track cut on distance from BarrelTracker inner r to determine whether track can form PFO"};
  Gaudi::Property<float> m_trackStateTolerance{
      this, "TrackStateTolerance", 0.0f,
      "Distance of possible second track state in the ECal Endcap to the ECal barrel inner radius"};
  Gaudi::Property<std::string> m_trackingSystemName{
      this, "TrackSystemName", "", "Name of the track fitting system to be used (KalTest, DDKalTest, aidaTT, ...)"};
  // For Strip Splitting method and also for hybrid ECAL
  Gaudi::Property<int> m_stripSplittingOn{this, "StripSplittingOn", 0,
                                          "To use strip splitting algorithm, this should be true"};
  // For Strip Splitting method and also for hybrid ECAL
  Gaudi::Property<int> m_useEcalScLayers{
      this, "UseEcalScLayers", {0}, "To use scintillator layers ~ hybrid ECAL, this should be true"};
  Gaudi::Property<int> m_useEcalSiLayers{
      this, "UseEcalSiLayers", {0}, "To use scintillator layers ~ hybrid ECAL, this should be true"};
  // Parameters for hybrid ECAL
  // Energy to MIP for Si-layers and Sc-layers, respectively.
  // Si
  Gaudi::Property<float> m_eCalSiToMip{this, "ECalSiToMipCalibration", 1.0f,
                                       "The calibration from deposited Si-layer energy to mip"};
  // Sc
  Gaudi::Property<float> m_eCalScToMip{this, "ECalScToMipCalibration", 1.0f,
                                       "The calibration from deposited Sc-layer energy to mip"};
  // MipThreshold for Si-layers and Sc-layers, respectively.
  // Si
  Gaudi::Property<float> m_eCalSiMipThreshold{this, "ECalSiMipThreshold", 0.0f,
                                              "Threshold for creating calo hits in the Si-layers of ECAL, units mip"};
  // Sc
  Gaudi::Property<float> m_eCalScMipThreshold{this, "ECalScMipThreshold", 0.0f,
                                              "Threshold for creating calo hits in the Sc-layers of ECAL, units mip"};
  // EcalToEM for Si-layers and Sc-layers, respectively.
  // Si
  Gaudi::Property<float> m_eCalSiToEMGeV{this, "ECalSiToEMGeVCalibration", 1.0f,
                                         "The calibration from deposited Si-layer energy to EM energy"};
  // Sc
  Gaudi::Property<float> m_eCalScToEMGeV{this, "ECalScToEMGeVCalibration", 1.0f,
                                         "The calibration from deposited Sc-layer energy to EM energy"};
  // EcalToHad for Si-layers and Sc-layers of the endcaps, respectively.
  // Si
  Gaudi::Property<float> m_eCalSiToHadGeVEndCap{
      this, "ECalSiToHadGeVCalibrationEndCap", 1.0f,
      "The calibration from deposited Si-layer energy on the endcaps to hadronic energy"};
  // Sc
  Gaudi::Property<float> m_eCalScToHadGeVEndCap{
      this, "ECalScToHadGeVCalibrationEndCap", 1.0f,
      "The calibration from deposited Sc-layer energy on the endcaps to hadronic energy"};
  // EcalToHad for Si-layers and Sc-layers of the barrel, respectively.
  // Si
  Gaudi::Property<float> m_eCalSiToHadGeVBarrel{
      this, "ECalSiToHadGeVCalibrationBarrel", 1.0f,
      "The calibration from deposited Si-layer energy on the barrel to hadronic energy"};
  // Sc
  Gaudi::Property<float> m_eCalScToHadGeVBarrel{
      this, "ECalScToHadGeVCalibrationBarrel", 1.0f,
      "The calibration from deposited Sc-layer energy to the barrel hadronic energy"};
  // Hadronic energy non-linearity correction
  Gaudi::Property<std::vector<float>> m_inputEnergyCorrectionPoints{
      this, "InputEnergyCorrectionPoints", {}, "The input energy points for hadronic energy correction"};
  Gaudi::Property<std::vector<float>> m_outputEnergyCorrectionPoints{
      this, "OutputEnergyCorrectionPoints", {}, "The output energy points for hadronic energy correction"};
  // ECAL energy non-linearity correction
  Gaudi::Property<std::vector<float>> m_ecalInputEnergyCorrectionPoints{
      this, "ECALInputEnergyCorrectionPoints", {}, "The input energy points for electromagnetic energy correction"};
  Gaudi::Property<std::vector<float>> m_ecalOutputEnergyCorrectionPoints{
      this, "ECALOutputEnergyCorrectionPoints", {}, "The output energy points for electromagnetic energy correction"};
  // EXTRA PARAMETERS FROM NIKIFOROS m_caloEncodingString
  Gaudi::Property<std::string> m_trackCreatorName{this, "TrackCreatorName", "DDTrackCreatorCLIC",
                                                  "The name of the DDTrackCreator implementation"};
  // Encoding Strings
  Gaudi::Property<std::string> m_trackEncodingStringVariable{
      this, "TrackerEncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for the Tracking detectors"};
  Gaudi::Property<std::string> m_caloEncodingStringVariable{
      this, "CaloEncodingStringParameterName", "GlobalCalorimeterReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for the Calo detectors"};
  /// EXTRA parameter that would initialize everything for ALLEGRO detector if m_detectorName=="ALLEGRO"
  Gaudi::Property<std::string> m_detectorName{this, "DetectorName", "", "The name of the detector"};
  Gaudi::Property<std::vector<float>> m_eCalBarrelNormalVector{
      this,
      "ECalBarrelNormalVector",
      {0.0f, 0.0f, 1.0f},
      "Normal vector for the ECal barrel sensitive layers in local coordinates"};
  Gaudi::Property<std::vector<float>> m_hCalBarrelNormalVector{
      this,
      "HCalBarrelNormalVector",
      {0.0f, 0.0f, 1.0f},
      "Normal vector for the HCal barrel sensitive layers in local coordinates"};
  Gaudi::Property<std::vector<float>> m_muonBarrelNormalVector{
      this,
      "YokeBarrelNormalVector",
      {0.0f, 0.0f, 1.0f},
      "Normal vector for the muon barrel sensitive layers in local coordinates"};

  // Reuse the software compensation parameters default values
  lc_content::LCSoftwareCompensationParameters defaultSoftwareCompensationParameters;
  Gaudi::Property<std::vector<float>> m_softCompParameters{
      this, "SoftwareCompensationWeights", defaultSoftwareCompensationParameters.m_softCompParameters,
      "The 9 software compensation weights for Pandora energy correction"};
  Gaudi::Property<std::vector<float>> m_softCompEnergyDensityBins{
      this, "SoftwareCompensationEnergyDensityBins", defaultSoftwareCompensationParameters.m_softCompEnergyDensityBins,
      "The 10 software compensation density bins for Pandora energy correction"};
  Gaudi::Property<float> m_energyDensityFinalBin{
      this, "FinalEnergyDensityBin", defaultSoftwareCompensationParameters.m_energyDensityFinalBin,
      "The software compensation final energy density bin for Pandora energy correction"};
  Gaudi::Property<float> m_maxClusterEnergyToApplySoftComp{
      this, "MaxClusterEnergyToApplySoftComp", defaultSoftwareCompensationParameters.m_maxClusterEnergyToApplySoftComp,
      "The maximum hadronic energy to apply software compensation in Pandora energy correction"};
  Gaudi::Property<float> m_minCleanHitEnergy{
      this, "MinCleanHitEnergy", defaultSoftwareCompensationParameters.m_minCleanHitEnergy,
      "The minimum hit energy to apply ECAL correction in Pandora energy correction"};
  Gaudi::Property<float> m_minCleanHitEnergyFraction{
      this, "MinCleanHitEnergyFraction", defaultSoftwareCompensationParameters.m_minCleanHitEnergyFraction,
      "The minimum hit energy fraction to apply ECAL correction in Pandora energy correction"};
  Gaudi::Property<float> m_minCleanCorrectedHitEnergy{
      this, "MinCleanCorrectedHitEnergy", defaultSoftwareCompensationParameters.m_minCleanCorrectedHitEnergy,
      "The minimum correction to an ECAL hit in Pandora energy correction"};
};

//------------------------------------------------------------------------------------------------------------------------------------------

#endif // #ifndef DDPANDORAPFANEWALGORITHM_H
