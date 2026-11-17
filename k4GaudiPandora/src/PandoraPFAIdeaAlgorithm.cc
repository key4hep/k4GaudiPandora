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
#include "PandoraPFAIdeaAlgorithm.h"
#include "DDExternalClusteringAlgorithm.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"

#include "LCClustering/EcalSeededClusteringAlgorithm.h"
#include "LCContent.h"
#include "LCParticleId/ForwardPhotonIdAlgorithm.h"
#include "LCPfoConstruction/IdeaPfoCreationAlgorithm.h"
#include "LCPlugins/DualReadoutCorrection.h"
#include "LCUtility/IsolatedHitPreparationAlgorithm.h"
#include "MLInference/ClusterNeutralPidAlgorithm.h"
#include "MLInference/SatelliteAssignmentOnnxAlgorithm.h"

#include "DDBFieldPlugin.h"
#include "DDPandoraPFANewAlgorithm.h"
#include "GeometryCreatorIdea.h"

PandoraPFAIdeaAlgorithm::PandoraPFAIdeaAlgorithm(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(
          name, svcLoc,
          {
              KeyValue("inputTrackCollection", "TracksFromGenParticles"),
              KeyValues("inputCaloHitCollections", {}),
              KeyValues("inputClusterCollections", {}),
          },
          {KeyValue("outputClusterCollection", "PandoraClusters"), KeyValue("outputPfoCollection", "PandoraPfaIdea")}),
      m_pandora() {}

StatusCode PandoraPFAIdeaAlgorithm::initialize() {
  // The pandora API and getExtension report failures by throwing; initialize() owes Gaudi a
  // StatusCode, so nothing is rethrown.
  try {
    m_geoSvc = serviceLocator()->service("GeoSvc");
    if (!m_geoSvc) {
      error() << "Unable to retrieve the GeoSvc" << endmsg;
      return StatusCode::FAILURE;
    }

    if (finaliseSteeringParameters().isFailure())
      return StatusCode::FAILURE;

    m_geometryCreator = std::make_unique<GeometryCreatorIdea>(m_geometryCreatorSettings, m_pandora, this);
    m_caloHitCreator = std::make_unique<DualReadoutCaloHitCreator>(m_caloHitCreatorSettings, m_pandora, this);
    m_trackCreator = std::make_unique<TrackCreatorIdea>(m_trackCreatorSettings, m_pandora, this);

    // TrackClusterAssociation, IsolatedHitMerging and VisualMonitoring come from here; only the
    // IDEA-specific algorithms are registered individually below.
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterAlgorithms(m_pandora))

    // Pandora takes ownership of the plugin, but the pointer is kept here so that chi can be read
    // back from it once the settings xml has been parsed - see below.
    auto* pDualReadoutCorrection = new lc_content::DualReadoutCorrection;
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::RegisterEnergyCorrectionPlugin(m_pandora, "DualReadoutCorrection",
                                                                       pandora::EnergyCorrectionType::HADRONIC,
                                                                       pDualReadoutCorrection));

    // Magnetic field from the dd4hep field map: algorithms retrieve it via the plugin (position
    // dependent) instead of a hardcoded XML value.
    dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::SetBFieldPlugin(m_pandora, new DDBFieldPlugin(mainDetector)));

    // Register algorithms
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::RegisterAlgorithmFactory(m_pandora, "DDExternalClustering",
                                                                 new DDExternalClusteringAlgorithm::Factory));

    // Set external parameters for DDExternalClusteringAlgorithm
    // ExternalClusterHolder is owned by this algorithm
    // ExternalEventParameter is created by this algo and deleted by Pandora
    m_extEvtParam = new ExternalEventParameter();
    m_extClusterHolder = std::make_unique<ExternalClusterHolder>();
    m_extEvtParam->m_externalClusterHolder = m_extClusterHolder.get();

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::SetExternalParameters(m_pandora, "DDExternalClustering", m_extEvtParam))

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::RegisterAlgorithmFactory(m_pandora, "ClusterNeutralPid",
                                                                 new lc_content::ClusterNeutralPidAlgorithm::Factory));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::RegisterAlgorithmFactory(m_pandora, "ForwardPhotonId",
                                                                 new lc_content::ForwardPhotonIdAlgorithm::Factory));

    PANDORA_THROW_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=,
        PandoraApi::RegisterAlgorithmFactory(m_pandora, "IsolatedHitPreparation",
                                             new lc_content::IsolatedHitPreparationAlgorithm::Factory));

    PANDORA_THROW_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=,
        PandoraApi::RegisterAlgorithmFactory(m_pandora, "EcalSeededClustering",
                                             new lc_content::EcalSeededClusteringAlgorithm::Factory));

    PANDORA_THROW_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=,
        PandoraApi::RegisterAlgorithmFactory(m_pandora, "SatelliteAssignmentOnnx",
                                             new lc_content::SatelliteAssignmentOnnxAlgorithm::Factory));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::RegisterAlgorithmFactory(m_pandora, "CreatePfo",
                                                                 new lc_content::IdeaPfoCreationAlgorithm::Factory));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_geometryCreator->CreateGeometry())

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::ReadSettings(m_pandora, m_pandoraSettingsXmlFile))

    // chi is only known once pandora has parsed the settings xml, so the pfo creator is built here
    // rather than alongside the other creators above.  Taking chi from the very plugin instance
    // that will apply the correction keeps the cluster energy error from drifting out of sync with
    // the energy itself.
    m_pfoCreatorSettings.m_chiEcal = pDualReadoutCorrection->GetChiEcal();
    m_pfoCreatorSettings.m_chiHcal = pDualReadoutCorrection->GetChiHcal();
    m_pfoCreator = std::make_unique<PfoCreatorIdea>(m_pfoCreatorSettings, m_pandora, this);

    return StatusCode::SUCCESS;
  } catch (const pandora::StatusCodeException& statusCodeException) {
    // pandora::StatusCodeException does NOT derive from std::exception, so it
    // must be caught explicitly (otherwise it falls through to catch(...)).
    error() << "Failed to initialize PandoraPFAIdeaAlgorithm: " << statusCodeException.ToString() << endmsg;
  } catch (const std::exception& exception) {
    error() << "Failed to initialize PandoraPFAIdeaAlgorithm: " << exception.what() << endmsg;
  } catch (...) {
    error() << "Failed to initialize PandoraPFAIdeaAlgorithm: unrecognized exception" << endmsg;
  }

  return StatusCode::FAILURE;
}

const pandora::Pandora* PandoraPFAIdeaAlgorithm::GetPandora() const { return &m_pandora; }

std::tuple<edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection>
PandoraPFAIdeaAlgorithm::operator()(const edm4hep::TrackCollection& trackColl,
                                    const std::vector<const edm4hep::CalorimeterHitCollection*>& caloHitColls,
                                    const std::vector<const edm4hep::ClusterCollection*>& clusterColls) const {

  try {
    // Create output collections
    edm4hep::ClusterCollection outClusterColl;
    edm4hep::ReconstructedParticleCollection pfoColl;

    // track
    std::vector<edm4hep::Track> tracksVector;

    for (const auto& aTrk : trackColl)
      tracksVector.push_back(aTrk);

    // TODO: track creator is working in progress
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_trackCreator->CreateTracks(tracksVector));

    // calo hits - create vectors to ensure stable addresses
    std::vector<std::vector<edm4hep::CalorimeterHit>> caloHitVectors(caloHitColls.size());
    for (size_t i = 0; i < caloHitColls.size(); ++i) {
      caloHitVectors[i].reserve(caloHitColls[i]->size());

      for (const auto& hit : *caloHitColls[i])
        caloHitVectors[i].push_back(hit);
    }

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, m_caloHitCreator->createCaloHits(caloHitVectors));

    // host edm4hep clusters for the external clustering algorithm
    std::unique_ptr<std::vector<std::vector<edm4hep::Cluster>>> externalClustersPtr =
        std::make_unique<std::vector<std::vector<edm4hep::Cluster>>>();
    externalClustersPtr->reserve(clusterColls.size());

    // loop over the input cluster collections and fill the external clusters vector
    for (const auto* clusterCollection : clusterColls) {
      std::vector<edm4hep::Cluster> clusterValues;
      clusterValues.reserve(clusterCollection->size());

      for (const auto& cluster : *clusterCollection) {
        clusterValues.push_back(cluster);
      }

      externalClustersPtr->push_back(std::move(clusterValues));
    }

    m_extClusterHolder->setExternalClusters(externalClustersPtr.get());

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(m_pandora));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            m_pfoCreator->CreateParticleFlowObjects(outClusterColl, pfoColl))

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(m_pandora))

    return std::make_tuple(std::move(outClusterColl), std::move(pfoColl));
  } catch (const pandora::StatusCodeException& statusCodeException) {
    // pandora::StatusCodeException does NOT derive from std::exception, so it
    // must be caught explicitly (otherwise it falls through to catch(...)).
    error() << "Pandora failed to process event: pandora::StatusCodeException " << statusCodeException.ToString()
            << endmsg;
    error() << statusCodeException.GetBackTrace() << endmsg;
    throw;
  } catch (std::exception& e) {
    error() << "Pandora failed to process event: std::exception " << e.what() << endmsg;
    throw;
  } catch (...) {
    error() << "Pandora failed to process event: unrecognized exception" << endmsg;
    throw;
  }
}

StatusCode PandoraPFAIdeaAlgorithm::finaliseSteeringParameters() {
  m_geometryCreatorSettings.m_hasHcalEndcap = m_hasHcalEndcap;
  if (!m_hasHcalEndcap) {
    warning() << "HasHcalEndcap is FALSE: no dual-readout HCAL endcap will be registered with "
                 "pandora, and every direction is treated as barrel.  This is only correct for a "
                 "barrel-only geometry for CI -- do NOT use it for production."
              << endmsg;
  }

  // TODO avoid duplication with DDPandoraPFANewAlgorithm
  auto getFieldFromCompact = []() -> double {
    dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
    const double position[3] = {0, 0, 0};      // position to calculate magnetic field at (the origin in this case)
    double magneticFieldVector[3] = {0, 0, 0}; // initialise object to hold magnetic field
    mainDetector.field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep

    return magneticFieldVector[2] / dd4hep::tesla; // z component at (0,0,0)
  };

  // The same selection GeometryCreatorIdea passes to SetEcalParameters, i.e. this is the SCEPCal
  // (ECAL) extension -- not the dual-readout one, despite the name it used to carry here.
  const dd4hep::rec::LayeredCalorimeterData* ecalExtension =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ENDCAP),
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

  // track creator settings
  m_trackCreatorSettings.m_bField = getFieldFromCompact();
  m_trackCreatorSettings.m_eCalEndCapInnerZ = ecalExtension->extent[2] / dd4hep::mm;
  // Needed by DDTrackCreatorBase::CalculateTrackTimeAtCalorimeter.  The SCEPCal barrel is
  // cylindrical, so m_eCalBarrelInnerSymmetry is left at 0 and only the radius is required.
  m_trackCreatorSettings.m_eCalBarrelInnerR = ecalExtension->extent[0] / dd4hep::mm;

  const size_t nSubDetectors = m_systemIDs.value().size();
  if (nSubDetectors == 0) {
    error() << "CaloSystemIDs is empty, no calorimeter is configured" << endmsg;
    return StatusCode::FAILURE;
  }

  auto checkLength = [this, nSubDetectors](const auto& prop) {
    if (prop.value().size() == nSubDetectors)
      return true;
    error() << prop.name() << " has " << prop.value().size() << " entries but CaloSystemIDs has " << nSubDetectors
            << "; the calo steering properties are indexed together, one entry per calorimeter" << endmsg;
    return false;
  };

  // single & rather than &&, so that every mismatched property is reported instead of just the first
  if (!(checkLength(m_collectionTypes) & checkLength(m_layerFieldNames) & checkLength(m_encodingStrings) &
        checkLength(m_cellSizes)))
    return StatusCode::FAILURE;

  // calo hit creator settings
  m_caloHitCreatorSettings.m_cherenkovFieldName = m_cherenkovFieldName;
  m_caloHitCreatorSettings.m_theta = std::atan2(ecalExtension->extent[0], ecalExtension->extent[2]);
  m_caloHitCreatorSettings.m_subDetectorSettings.resize(nSubDetectors);

  for (size_t icol = 0; icol < nSubDetectors; ++icol) {
    auto& subdetectorSetting = m_caloHitCreatorSettings.m_subDetectorSettings.at(icol);
    subdetectorSetting.m_systemID = m_systemIDs.value().at(icol);
    subdetectorSetting.m_encodingString = m_encodingStrings.value().at(icol);
    subdetectorSetting.m_layerFieldName = m_layerFieldNames.value().at(icol);
    subdetectorSetting.m_collectionType = m_collectionTypes.value().at(icol);
    subdetectorSetting.m_cellSize = m_cellSizes.value().at(icol);

    // A dual-readout tube is a single channel spanning the full depth, so a subdetector with no
    // layer field has no per-layer thickness and DualReadoutCaloHitCreator falls back to the cell size.
    if (!subdetectorSetting.m_layerFieldName.empty())
      for (const auto& layer : ecalExtension->layers)
        subdetectorSetting.m_layerThicknesses.push_back(layer.sensitive_thickness);
  }

  return StatusCode::SUCCESS;
}

DECLARE_COMPONENT(PandoraPFAIdeaAlgorithm)
