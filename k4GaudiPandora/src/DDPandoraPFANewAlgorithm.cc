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
 *  @file   DDMarlinPandora/src/DDPandoraPFANewAlgorithm.cc
 *
 *  @brief  Implementation of the pandora pfa new processor class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "LCContent.h"
#include "LCPlugins/LCSoftwareCompensation.h"


#include "DDPandoraPFANewAlgorithm.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DDRec/DetectorData.h"

#include "DDTrackCreatorCLIC.h"
//#include "DDTrackCreatorILD.h"

#include "DDBFieldPlugin.h"

#include <cstdlib>

DECLARE_COMPONENT(DDPandoraPFANewAlgorithm)

double getFieldFromCompact() {
  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
  const double      position[3]  = {0, 0, 0};  // position to calculate magnetic field at (the origin in this case)
  double            magneticFieldVector[3] = {0, 0, 0};               // initialise object to hold magnetic field
  mainDetector.field().magneticField(position, magneticFieldVector);  // get the magnetic field vector from DD4hep

  return magneticFieldVector[2] / dd4hep::tesla;  // z component at (0,0,0)
}

dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, MsgStream log,unsigned int excludeFlag = 0) {
  dd4hep::rec::LayeredCalorimeterData* theExtension = 0;

  dd4hep::Detector&                      mainDetector = dd4hep::Detector::getInstance();
  const std::vector<dd4hep::DetElement>& theDetectors =
      dd4hep::DetectorSelector(mainDetector).detectors(includeFlag, excludeFlag);

  log << MSG::DEBUG << " getExtension :  includeFlag: " << dd4hep::DetType(includeFlag)
                    << " excludeFlag: " << dd4hep::DetType(excludeFlag) << "  found : " << theDetectors.size()
                    << "  - first det: " << theDetectors.at(0).name() << endmsg;

  if (theDetectors.size() != 1) {
    std::stringstream es;
    es << " getExtension: selection is not unique (or empty)  includeFlag: " << dd4hep::DetType(includeFlag)
       << " excludeFlag: " << dd4hep::DetType(excludeFlag) << " --- found detectors : ";
    for (unsigned i = 0, N = theDetectors.size(); i < N; ++i) {
      es << theDetectors.at(i).name() << ", ";
    }
    throw std::runtime_error(es.str());
  }

  theExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

  return theExtension;
}

std::vector<double> getTrackingRegionExtent() {
  ///Rmin, Rmax, Zmax
  std::vector<double> extent;

  extent.reserve(3);

  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();

  extent[0] = 0.1;  ///FIXME! CLIC-specific: Inner radius was set to 0 for SiD-type detectors
  extent[1] = mainDetector.constantAsDouble("tracker_region_rmax") / dd4hep::mm;
  extent[2] = mainDetector.constantAsDouble("tracker_region_zmax") / dd4hep::mm;

  return extent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDPandoraPFANewAlgorithm::DDPandoraPFANewAlgorithm(const std::string& name, ISvcLocator* svcLoc) : MultiTransformer(name, svcLoc,
  { KeyValues("MCParticleCollections", {"MCParticles"}), 
    KeyValues("KinkVertexCollections", {}),
    KeyValues("ProngVertexCollections", {}),
    KeyValues("SplitVertexCollections", {}),
    KeyValues("V0VertexCollections", {}),
    KeyValues("TrackCollections", {"TrackCollections"}),    
    KeyValues("RelTrackCollections", {"RelTrackCollections"}),
    KeyValues("ECalCaloHitCollections", {"ECalCaloHitCollections"}),
    KeyValues("HCalCaloHitCollections", {"HCalCaloHitCollections"}),
    KeyValues("MuonCaloHitCollections", {}),
    KeyValues("LCalCaloHitCollections", {}),
    KeyValues("LHCalCaloHitCollections", {}),
    KeyValues("RelCaloHitCollections", {"RelCaloHitCollections"}),},
  { KeyValues("ClusterCollectionName", {"PandoraPFANewClusters"}),
    KeyValues("PFOCollectionName", {"PandoraPFANewPFOs"}), 
    KeyValues("StartVertexCollectionName", {"PandoraPFANewStartVertices"}) }) {}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DDPandoraPFANewAlgorithm::initialize() {
  m_geoSvc = serviceLocator()->service("GeoSvc");  // important to initialize m_geoSvc
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  try {
    info() << "DDPandoraPFANewAlgorithm - Init" << endmsg;
    this->FinaliseSteeringParameters();

    m_pPandora         = new pandora::Pandora();
    m_pGeometryCreator = new DDGeometryCreator(m_geometryCreatorSettings, m_pPandora, msgSvc());
    m_pCaloHitCreator  = new DDCaloHitCreator(m_caloHitCreatorSettings, m_pPandora, msgSvc());

    ///FIXME: IMPLEMENT FACTORY
    if (m_settings.m_trackCreatorName == "DDTrackCreatorCLIC")
      m_pTrackCreator = new DDTrackCreatorCLIC(m_trackCreatorSettings, m_pPandora, msgSvc());
    //else if (m_settings.m_trackCreatorName == "DDTrackCreatorILD")
      //m_pTrackCreator = new DDTrackCreatorILD(m_trackCreatorSettings, m_pPandora);
    else
      error() << "Unknown DDTrackCreator: " << m_settings.m_trackCreatorName << endmsg;

    m_pDDMCParticleCreator = new DDMCParticleCreator(m_mcParticleCreatorSettings, m_pPandora, msgSvc());
    m_pDDPfoCreator        = new DDPfoCreator(m_pfoCreatorSettings, m_pPandora, msgSvc());

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->RegisterUserComponents());
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pGeometryCreator->CreateGeometry());
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::ReadSettings(*m_pPandora, m_settings.m_pandoraSettingsXmlFile));

  } catch (pandora::StatusCodeException& statusCodeException) {
    error() << "Failed to initialize marlin pandora: " << statusCodeException.ToString() << endmsg;
    throw statusCodeException;
  } catch (std::exception& exception) {
    error() << "Failed to initialize marlin pandora: std exception " << exception.what() << endmsg;
    throw exception;
  } catch (...) {
    error() << "Failed to initialize marlin pandora: unrecognized exception" << endmsg;
    throw;
  }

  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::tuple<edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection, edm4hep::VertexCollection> DDPandoraPFANewAlgorithm::operator()(
  const std::vector<const edm4hep::MCParticleCollection*>& MCParticleCollections,
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
  const std::vector<const edm4hep::CaloHitSimCaloHitLinkCollection*>& caloLinkCollections
) const {
  try {
    debug() << "DDPandoraPFANewAlgorithm - Run " << endmsg;
    //(void)m_pandoraToLCEventMap.insert(PandoraToLCEventMap::value_type(m_pPandora, pLCEvent));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pDDMCParticleCreator->CreateMCParticles(MCParticleCollections));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTrackAssociations(kinkCollections, prongCollections, splitCollections, v0Collections));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_pTrackCreator->CreateTracks(trackCollections));
    PANDORA_THROW_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=,
        m_pDDMCParticleCreator->CreateTrackToMCParticleRelationships(trackerHitLinkCollections, m_pTrackCreator->GetTrackVector()));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, 
      m_pCaloHitCreator->CreateCaloHits(eCalCollections, hCalCollections, mCalCollections, lCalCollections, lhCalCollections));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            m_pDDMCParticleCreator->CreateCaloHitToMCParticleRelationships(
                              caloLinkCollections, m_pCaloHitCreator->GetCalorimeterHitVector()));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));

    edm4hep::ClusterCollection               pClusterCollection;
    edm4hep::ReconstructedParticleCollection pReconstructedParticleCollection;
    edm4hep::VertexCollection                pStartVertexCollection;

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, 
      m_pDDPfoCreator->CreateParticleFlowObjects(pClusterCollection, pReconstructedParticleCollection, pStartVertexCollection));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
    this->Reset();

    return std::make_tuple(
      std::move(pClusterCollection), 
      std::move(pReconstructedParticleCollection), 
      std::move(pStartVertexCollection)
    );
  } catch (pandora::StatusCodeException& statusCodeException) {
    error() << "Marlin pandora failed to process event: " << statusCodeException.ToString() << endmsg;
    throw statusCodeException;
  } catch (std::exception& exception) {
    error() << "Pandora failed to process event: std exception " << exception.what() << endmsg;
    throw;
  } catch (...) {
    error() << "Pandora failed to process event: unrecognized exception" << endmsg;
    throw;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DDPandoraPFANewAlgorithm::finalize() {
  delete m_pPandora;
  delete m_pGeometryCreator;
  delete m_pCaloHitCreator;
  delete m_pTrackCreator;
  delete m_pDDMCParticleCreator;
  delete m_pDDPfoCreator;

  info() << "DDPandoraPFANewAlgorithm - End" << endmsg;
  return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora* DDPandoraPFANewAlgorithm::GetPandora() const {
  if (NULL == m_pPandora)
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

  return m_pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------
/*
const EVENT::LCEvent* DDPandoraPFANewAlgorithm::GetCurrentEvent(const pandora::Pandora* const pPandora) {
  PandoraToLCEventMap::iterator iter = m_pandoraToLCEventMap.find(pPandora);

  if (m_pandoraToLCEventMap.end() == iter)
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

  return iter->second;
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDPandoraPFANewAlgorithm::RegisterUserComponents() const {
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterAlgorithms(*m_pPandora));
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterBasicPlugins(*m_pPandora));

  if (m_settings.m_useDD4hepField) {
    dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                             PandoraApi::SetBFieldPlugin(*m_pPandora, new DDBFieldPlugin(mainDetector)));
  } else {
    PANDORA_RETURN_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=,
        LCContent::RegisterBFieldPlugin(*m_pPandora, m_settings.m_innerBField, m_settings.m_muonBarrelBField,
                                        m_settings.m_muonEndCapBField));
  }

  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                           LCContent::RegisterNonLinearityEnergyCorrection(
                               *m_pPandora, "NonLinearity", pandora::HADRONIC, m_settings.m_inputEnergyCorrectionPoints,
                               m_settings.m_outputEnergyCorrectionPoints));

  lc_content::LCSoftwareCompensationParameters softwareCompensationParameters;
  softwareCompensationParameters.m_softCompParameters              = m_settings.m_softCompParameters;
  softwareCompensationParameters.m_softCompEnergyDensityBins       = m_settings.m_softCompEnergyDensityBins;
  softwareCompensationParameters.m_energyDensityFinalBin           = m_settings.m_energyDensityFinalBin;
  softwareCompensationParameters.m_maxClusterEnergyToApplySoftComp = m_settings.m_maxClusterEnergyToApplySoftComp;
  softwareCompensationParameters.m_minCleanHitEnergy               = m_settings.m_minCleanHitEnergy;
  softwareCompensationParameters.m_minCleanHitEnergyFraction       = m_settings.m_minCleanHitEnergyFraction;
  softwareCompensationParameters.m_minCleanCorrectedHitEnergy      = m_settings.m_minCleanCorrectedHitEnergy;

  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                           LCContent::RegisterSoftwareCompensationEnergyCorrection(*m_pPandora, "SoftwareCompensation",
                                                                                   softwareCompensationParameters));

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewAlgorithm::FinaliseSteeringParameters() {
  // ATTN: This function seems to be necessary for operations that cannot easily be performed at construction of the processor,
  // when the steering file is parsed e.g. the call to GEAR to get the inner bfield
  MsgStream log(msgSvc(), name());

  m_settings.m_pandoraSettingsXmlFile = m_pandoraSettingsXmlFile;
  m_geometryCreatorSettings.m_createGaps = m_createGaps;
  m_pfoCreatorSettings.m_startVertexAlgName = m_startVertexAlgName;
  m_pfoCreatorSettings.m_emStochasticTerm = m_emStochasticTerm;
  m_pfoCreatorSettings.m_hadStochasticTerm = m_hadStochasticTerm;
  m_pfoCreatorSettings.m_emConstantTerm = m_emConstantTerm;
  m_pfoCreatorSettings.m_hadConstantTerm = m_hadConstantTerm;
  m_caloHitCreatorSettings.m_eCalToMip = m_eCalToMip;
  m_caloHitCreatorSettings.m_hCalToMip = m_hCalToMip;
  m_caloHitCreatorSettings.m_eCalMipThreshold = m_eCalMipThreshold;
  m_caloHitCreatorSettings.m_muonToMip = m_muonToMip;
  m_caloHitCreatorSettings.m_hCalMipThreshold = m_hCalMipThreshold;
  m_caloHitCreatorSettings.m_eCalToEMGeV = m_eCalToEMGeV;
  m_caloHitCreatorSettings.m_hCalToEMGeV = m_hCalToEMGeV;
  m_caloHitCreatorSettings.m_eCalToHadGeVEndCap = m_eCalToHadGeVEndCap;
  m_caloHitCreatorSettings.m_eCalToHadGeVBarrel = m_eCalToHadGeVBarrel;
  m_caloHitCreatorSettings.m_hCalToHadGeV = m_hCalToHadGeV;
  m_caloHitCreatorSettings.m_muonDigitalHits = m_muonDigitalHits;
  m_caloHitCreatorSettings.m_muonHitEnergy = m_muonHitEnergy;
  m_caloHitCreatorSettings.m_maxHCalHitHadronicEnergy = m_maxHCalHitHadronicEnergy;
  m_caloHitCreatorSettings.m_nOuterSamplingLayers = m_nOuterSamplingLayers;
  m_caloHitCreatorSettings.m_layersFromEdgeMaxRearDistance = m_layersFromEdgeMaxRearDistance;
  m_settings.m_muonBarrelBField = m_muonBarrelBField;
  m_settings.m_muonEndCapBField = m_muonEndCapBField;
  m_settings.m_useDD4hepField = m_useDD4hepField;
  m_trackCreatorSettings.m_shouldFormTrackRelationships = m_shouldFormTrackRelationships;
  m_trackCreatorSettings.m_minTrackHits = m_minTrackHits;
  m_trackCreatorSettings.m_minFtdTrackHits = m_minFtdTrackHits;
  m_trackCreatorSettings.m_maxTrackHits = m_maxTrackHits;
  m_trackCreatorSettings.m_d0TrackCut = m_d0TrackCut;
  m_trackCreatorSettings.m_z0TrackCut = m_z0TrackCut;
  m_trackCreatorSettings.m_usingNonVertexTracks = m_usingNonVertexTracks;
  m_trackCreatorSettings.m_usingUnmatchedNonVertexTracks = m_usingUnmatchedNonVertexTracks;
  m_trackCreatorSettings.m_usingUnmatchedVertexTracks = m_usingUnmatchedVertexTracks;
  m_trackCreatorSettings.m_unmatchedVertexTrackMaxEnergy = m_unmatchedVertexTrackMaxEnergy;
  m_trackCreatorSettings.m_d0UnmatchedVertexTrackCut = m_d0UnmatchedVertexTrackCut;
  m_trackCreatorSettings.m_z0UnmatchedVertexTrackCut = m_z0UnmatchedVertexTrackCut;
  m_trackCreatorSettings.m_zCutForNonVertexTracks = m_zCutForNonVertexTracks;
  m_trackCreatorSettings.m_reachesECalNBarrelTrackerHits = m_reachesECalNBarrelTrackerHits;
  m_trackCreatorSettings.m_reachesECalNFtdHits = m_reachesECalNFtdHits;
  m_trackCreatorSettings.m_reachesECalBarrelTrackerOuterDistance = m_reachesECalBarrelTrackerOuterDistance;
  m_trackCreatorSettings.m_reachesECalMinFtdLayer = m_reachesECalMinFtdLayer;
  m_trackCreatorSettings.m_reachesECalBarrelTrackerZMaxDistance = m_reachesECalBarrelTrackerZMaxDistance;
  m_trackCreatorSettings.m_reachesECalFtdZMaxDistance = m_reachesECalFtdZMaxDistance;
  m_trackCreatorSettings.m_curvatureToMomentumFactor = m_curvatureToMomentumFactor;
  m_trackCreatorSettings.m_minTrackECalDistanceFromIp = m_minTrackECalDistanceFromIp;
  m_trackCreatorSettings.m_maxTrackSigmaPOverP = m_maxTrackSigmaPOverP;
  m_trackCreatorSettings.m_minMomentumForTrackHitChecks = m_minMomentumForTrackHitChecks;
  m_trackCreatorSettings.m_minBarrelTrackerHitFractionOfExpected = m_minBarrelTrackerHitFractionOfExpected;
  m_trackCreatorSettings.m_minFtdHitsForBarrelTrackerHitFraction = m_minFtdHitsForBarrelTrackerHitFraction;
  m_trackCreatorSettings.m_maxBarrelTrackerInnerRDistance = m_maxBarrelTrackerInnerRDistance;
  m_trackCreatorSettings.m_trackStateTolerance = m_trackStateTolerance;
  m_trackCreatorSettings.m_trackingSystemName = m_trackingSystemName;
  std::string trackString;
  trackString = m_geoSvc->constantAsString(m_trackEncodingStringVariable.value());
  m_trackCreatorSettings.m_trackingEncodingString = trackString;
  m_caloHitCreatorSettings.m_stripSplittingOn = m_stripSplittingOn;
  m_caloHitCreatorSettings.m_useEcalScLayers = m_useEcalScLayers;
  m_caloHitCreatorSettings.m_useEcalSiLayers = m_useEcalSiLayers;
  m_caloHitCreatorSettings.m_eCalSiToMip = m_eCalSiToMip;
  m_caloHitCreatorSettings.m_eCalScToMip = m_eCalScToMip;
  m_caloHitCreatorSettings.m_eCalSiMipThreshold = m_eCalSiMipThreshold;
  m_caloHitCreatorSettings.m_eCalScMipThreshold = m_eCalScMipThreshold;
  m_caloHitCreatorSettings.m_eCalSiToEMGeV = m_eCalSiToEMGeV;
  m_caloHitCreatorSettings.m_eCalScToEMGeV = m_eCalScToEMGeV;
  m_caloHitCreatorSettings.m_eCalSiToHadGeVEndCap = m_eCalSiToHadGeVEndCap;
  m_caloHitCreatorSettings.m_eCalScToHadGeVEndCap = m_eCalScToHadGeVEndCap;
  m_caloHitCreatorSettings.m_eCalSiToHadGeVBarrel = m_eCalSiToHadGeVBarrel;
  m_caloHitCreatorSettings.m_eCalScToHadGeVBarrel = m_eCalScToHadGeVBarrel;
  std::string caloString;
  caloString = m_geoSvc->constantAsString(m_caloEncodingStringVariable.value());
  m_caloHitCreatorSettings.m_caloEncodingString = caloString;
  m_settings.m_inputEnergyCorrectionPoints = m_inputEnergyCorrectionPoints;
  m_settings.m_outputEnergyCorrectionPoints = m_outputEnergyCorrectionPoints;
  m_settings.m_ecalInputEnergyCorrectionPoints = m_ecalInputEnergyCorrectionPoints;
  m_settings.m_ecalOutputEnergyCorrectionPoints = m_ecalOutputEnergyCorrectionPoints;
  m_settings.m_trackCreatorName = m_trackCreatorName;
  m_settings.m_detectorName = m_detectorName;
  m_caloHitCreatorSettings.m_eCalBarrelNormalVector = m_eCalBarrelNormalVector;
  m_caloHitCreatorSettings.m_hCalBarrelNormalVector = m_hCalBarrelNormalVector;
  m_caloHitCreatorSettings.m_muonBarrelNormalVector = m_muonBarrelNormalVector;
  m_settings.m_softCompParameters = m_softCompParameters;
  m_settings.m_softCompEnergyDensityBins = m_softCompEnergyDensityBins;
  m_settings.m_energyDensityFinalBin = m_energyDensityFinalBin;
  m_settings.m_maxClusterEnergyToApplySoftComp = m_maxClusterEnergyToApplySoftComp;
  m_settings.m_minCleanHitEnergy = m_minCleanHitEnergy;
  m_settings.m_minCleanHitEnergyFraction = m_minCleanHitEnergyFraction;
  m_settings.m_minCleanCorrectedHitEnergy = m_minCleanCorrectedHitEnergy;

  m_trackCreatorSettings.m_bField = getFieldFromCompact();

  //Get ECal Barrel extension by type, ignore plugs and rings
  const dd4hep::rec::LayeredCalorimeterData* eCalBarrelExtension =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL), log,
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));
  //Get ECal Endcap extension by type, ignore plugs and rings
  const dd4hep::rec::LayeredCalorimeterData* eCalEndcapExtension =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP), log,
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));
  //Get HCal Barrel extension by type, ignore plugs and rings
  const dd4hep::rec::LayeredCalorimeterData* hCalBarrelExtension =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL), log,
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));
  //Get HCal Endcap extension by type, ignore plugs and rings
  const dd4hep::rec::LayeredCalorimeterData* hCalEndcapExtension =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP), log,
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));
  //Get Muon Barrel extension by type, ignore plugs and rings
  const dd4hep::rec::LayeredCalorimeterData* muonBarrelExtension =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::BARREL), log,
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));
  //fg: muon endcap is not used :
  // //Get Muon Endcap extension by type, ignore plugs and rings
  // const dd4hep::rec::LayeredCalorimeterData * muonEndcapExtension= getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::ENDCAP),  log, ( dd4hep::DetType::AUXILIARY ) );

  //Get COIL extension
  const dd4hep::rec::LayeredCalorimeterData* coilExtension = getExtension((dd4hep::DetType::COIL), log);

  m_trackCreatorSettings.m_eCalBarrelInnerSymmetry = eCalBarrelExtension->inner_symmetry;
  m_trackCreatorSettings.m_eCalBarrelInnerPhi0     = eCalBarrelExtension->inner_phi0 / dd4hep::rad;
  m_trackCreatorSettings.m_eCalBarrelInnerR        = eCalBarrelExtension->extent[0] / dd4hep::mm;
  m_trackCreatorSettings.m_eCalEndCapInnerZ        = eCalEndcapExtension->extent[2] / dd4hep::mm;

  m_caloHitCreatorSettings.m_eCalBarrelOuterZ             = eCalBarrelExtension->extent[3] / dd4hep::mm;
  m_caloHitCreatorSettings.m_hCalBarrelOuterZ             = hCalBarrelExtension->extent[3] / dd4hep::mm;
  m_caloHitCreatorSettings.m_muonBarrelOuterZ             = muonBarrelExtension->extent[3] / dd4hep::mm;
  m_caloHitCreatorSettings.m_coilOuterR                   = coilExtension->extent[1] / dd4hep::mm;
  m_caloHitCreatorSettings.m_eCalBarrelInnerPhi0          = eCalBarrelExtension->inner_phi0 / dd4hep::rad;
  m_caloHitCreatorSettings.m_eCalBarrelInnerSymmetry      = eCalBarrelExtension->inner_symmetry;
  m_caloHitCreatorSettings.m_hCalBarrelInnerPhi0          = hCalBarrelExtension->inner_phi0 / dd4hep::rad;
  m_caloHitCreatorSettings.m_hCalBarrelInnerSymmetry      = hCalBarrelExtension->inner_symmetry;
  m_caloHitCreatorSettings.m_muonBarrelInnerPhi0          = muonBarrelExtension->inner_phi0 / dd4hep::rad;
  m_caloHitCreatorSettings.m_muonBarrelInnerSymmetry      = muonBarrelExtension->inner_symmetry;
  m_caloHitCreatorSettings.m_hCalEndCapOuterR             = hCalEndcapExtension->extent[1] / dd4hep::mm;
  m_caloHitCreatorSettings.m_hCalEndCapOuterZ             = hCalEndcapExtension->extent[3] / dd4hep::mm;
  m_caloHitCreatorSettings.m_hCalBarrelOuterR             = hCalBarrelExtension->extent[1] / dd4hep::mm;
  m_caloHitCreatorSettings.m_hCalBarrelOuterPhi0          = hCalBarrelExtension->outer_phi0 / dd4hep::rad;
  m_caloHitCreatorSettings.m_hCalBarrelOuterSymmetry      = hCalBarrelExtension->outer_symmetry;
  m_caloHitCreatorSettings.m_hCalEndCapInnerSymmetryOrder = hCalEndcapExtension->inner_symmetry;
  ;
  m_caloHitCreatorSettings.m_hCalEndCapInnerPhiCoordinate = hCalEndcapExtension->inner_phi0 / dd4hep::rad;
  ;

  // Get the magnetic field
  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
  const double      position[3]  = {0, 0, 0};  // position to calculate magnetic field at (the origin in this case)
  double            magneticFieldVector[3] = {0, 0, 0};               // initialise object to hold magnetic field
  mainDetector.field().magneticField(position, magneticFieldVector);  // get the magnetic field vector from DD4hep

  m_settings.m_innerBField = magneticFieldVector[2] / dd4hep::tesla;  // z component at (0,0,0)
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDPandoraPFANewAlgorithm::Reset() const{
  m_pCaloHitCreator->Reset();
  m_pTrackCreator->Reset();
  /*
  PandoraToLCEventMap::iterator iter = m_pandoraToLCEventMap.find(m_pPandora);

  if (m_pandoraToLCEventMap.end() == iter)
    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

  m_pandoraToLCEventMap.erase(iter);*/
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDPandoraPFANewAlgorithm::Settings::Settings()
    : m_innerBField(3.5f),
      m_muonBarrelBField(-1.5f),
      m_muonEndCapBField(0.01f),
      m_inputEnergyCorrectionPoints(0),
      m_outputEnergyCorrectionPoints(0),
      m_trackCreatorName("")

{}
