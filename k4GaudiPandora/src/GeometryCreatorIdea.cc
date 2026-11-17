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
 *  @file   k4GaudiPandora/src/GeometryCreatorIdea.cc
 *
 *  @brief  Implementation of the geometry creator class.
 *
 *  $Log: $
 */

#include "GeometryCreatorIdea.h"

#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"

// Forward declarations. See DDPandoraPFANewAlgorithm.cc
// dd4hep::rec::LayeredCalorimeterData * getExtension(std::string detectorName);
dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag = 0);

GeometryCreatorIdea::GeometryCreatorIdea(const Settings& settings, pandora::Pandora& pPandora,
                                         Gaudi::Algorithm* algorithm)
    : DDGeometryCreator(settings, pPandora, algorithm), m_hasHcalEndcap(settings.m_hasHcalEndcap) {}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode GeometryCreatorIdea::CreateGeometry() const {
  try {
    SubDetectorTypeMap subDetectorTypeMap;
    this->SetMandatorySubDetectorParameters(subDetectorTypeMap);

    m_algorithm.debug() << "Creating geometry for IDEA detector" << endmsg;

    for (SubDetectorTypeMap::const_iterator iter = subDetectorTypeMap.begin(), iterEnd = subDetectorTypeMap.end();
         iter != iterEnd; ++iter) {
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                               PandoraApi::Geometry::SubDetector::Create(m_pPandora, iter->second));
    }
  } catch (std::exception& exception) {
    m_algorithm.error() << "Failure in GeometryCreatorIdea, exception: " << exception.what() << endmsg;
    throw exception;
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GeometryCreatorIdea::SetMandatorySubDetectorParameters(SubDetectorTypeMap& subDetectorTypeMap) const {
  PandoraApi::Geometry::SubDetector::Parameters eCalBarrelParameters, eCalEndCapParameters, hCalBarrelParameters,
      hCalEndCapParameters;
  // hCalBarrelParameters, hCalEndCapParameters, muonBarrelParameters, muonEndCapParameters;
  // TODO they're not used anywhere at the moment, so ignoring them

  this->SetEcalParameters(*const_cast<dd4hep::rec::LayeredCalorimeterData*>(getExtension(
                              (dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ENDCAP),
                              (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))),
                          eCalBarrelParameters, eCalEndCapParameters);

  subDetectorTypeMap[pandora::ECAL_BARREL] = eCalBarrelParameters;
  subDetectorTypeMap[pandora::ECAL_ENDCAP] = eCalEndCapParameters;

  this->SetHcalBarrelParameters(
      *const_cast<dd4hep::rec::LayeredCalorimeterData*>(
          getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::HADRONIC),
                       (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))),
      hCalBarrelParameters);

  subDetectorTypeMap[pandora::HCAL_BARREL] = hCalBarrelParameters;

  // A barrel-only geometry has no dual-readout endcap to look up, and registering one would fail.
  if (m_hasHcalEndcap) {
    this->SetHcalEndcapParameters(
        *const_cast<dd4hep::rec::LayeredCalorimeterData*>(
            getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ENDCAP | dd4hep::DetType::HADRONIC),
                         (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))),
        hCalEndCapParameters);

    subDetectorTypeMap[pandora::HCAL_ENDCAP] = hCalEndCapParameters;
  }

  // PandoraApi::Geometry::SubDetector::Parameters coilParameters; // TODO retrive coil parameters
}

void GeometryCreatorIdea::SetEcalParameters(const dd4hep::rec::LayeredCalorimeterData& inputParameters,
                                            PandoraApi::Geometry::SubDetector::Parameters& paramBarrel,
                                            PandoraApi::Geometry::SubDetector::Parameters& paramEndcap) const {
  unsigned layerVecSize = inputParameters.layers.size();  // zero for option 1
  unsigned nlayers = layerVecSize > 0 ? layerVecSize : 1; // avoid zero

  paramBarrel.m_subDetectorName = "ECalBarrel";
  paramBarrel.m_subDetectorType = pandora::ECAL_BARREL;
  paramBarrel.m_innerRCoordinate = inputParameters.extent[0] / dd4hep::mm;
  paramBarrel.m_innerZCoordinate = 0.;
  paramBarrel.m_innerPhiCoordinate = 0.; // not initialized in the LayeredCalorimeterData
  paramBarrel.m_innerSymmetryOrder = 0;  // not initialized
  paramBarrel.m_outerRCoordinate =
      (inputParameters.extent[0] + inputParameters.extent[3] - inputParameters.extent[2]) / dd4hep::mm;
  // barrel outerR = barrel innerR + tower height, tower height = endcap outer Z - endcap inner Z
  paramBarrel.m_outerZCoordinate = inputParameters.extent[2] / dd4hep::mm; // use endcap inner Z
  paramBarrel.m_outerPhiCoordinate = 0.;                                   // not initialized
  paramBarrel.m_outerSymmetryOrder = 0;                                    // not initialized
  paramBarrel.m_isMirroredInZ = true;
  paramBarrel.m_nLayers = nlayers; // no longitudinal segmentation

  // just dummy values for the mandatory parameters
  paramBarrel.m_layerParametersVector.resize(nlayers);
  float distanceBarrel = inputParameters.extent[0];

  for (unsigned iLayer = 0; iLayer < nlayers; ++iLayer) {
    paramBarrel.m_layerParametersVector.at(iLayer).m_closestDistanceToIp = distanceBarrel / dd4hep::mm;
    paramBarrel.m_layerParametersVector.at(iLayer).m_nRadiationLengths = 0.;   // not used
    paramBarrel.m_layerParametersVector.at(iLayer).m_nInteractionLengths = 0.; // not used

    if (layerVecSize > 0)
      distanceBarrel += inputParameters.layers.at(iLayer).sensitive_thickness;
  }

  paramEndcap.m_subDetectorName = "ECalEndcap";
  paramEndcap.m_subDetectorType = pandora::ECAL_ENDCAP;
  paramEndcap.m_innerRCoordinate = inputParameters.extent[4] / dd4hep::mm;
  paramEndcap.m_innerZCoordinate = inputParameters.extent[2] / dd4hep::mm;
  paramEndcap.m_innerPhiCoordinate = 0.; // not initialized in the LayeredCalorimeterData
  paramEndcap.m_innerSymmetryOrder = 0;  // not initialized
  paramEndcap.m_outerRCoordinate = inputParameters.extent[5] / dd4hep::mm;
  paramEndcap.m_outerZCoordinate = inputParameters.extent[3] / dd4hep::mm;
  paramEndcap.m_outerPhiCoordinate = 0.; // not initialized
  paramEndcap.m_outerSymmetryOrder = 0;  // not initialized
  paramEndcap.m_isMirroredInZ = true;
  paramEndcap.m_nLayers = nlayers; // no longitudinal segmentation

  // just dummy values for the mandatory parameters
  paramEndcap.m_layerParametersVector.resize(nlayers);
  float distanceEndcap = inputParameters.extent[2];

  for (unsigned iLayer = 0; iLayer < nlayers; ++iLayer) {
    paramEndcap.m_layerParametersVector.at(iLayer).m_closestDistanceToIp = distanceEndcap / dd4hep::mm;
    paramEndcap.m_layerParametersVector.at(iLayer).m_nRadiationLengths = 0.;   // not used
    paramEndcap.m_layerParametersVector.at(iLayer).m_nInteractionLengths = 0.; // not used

    if (layerVecSize > 0)
      distanceEndcap += inputParameters.layers.at(iLayer).sensitive_thickness;
  }

  return;
}

void GeometryCreatorIdea::SetHcalBarrelParameters(const dd4hep::rec::LayeredCalorimeterData& inputParameters,
                                                  PandoraApi::Geometry::SubDetector::Parameters& paramBarrel) const {

  unsigned layerVecSize = 0.;
  unsigned nlayers = layerVecSize > 0 ? layerVecSize : 1; // avoid zero

  paramBarrel.m_subDetectorName = "HCalBarrel";
  paramBarrel.m_subDetectorType = pandora::HCAL_BARREL;
  paramBarrel.m_innerRCoordinate = inputParameters.extent[0] / dd4hep::mm;
  paramBarrel.m_innerZCoordinate = 0.;
  paramBarrel.m_innerPhiCoordinate = 0.;
  paramBarrel.m_innerSymmetryOrder = 0;
  paramBarrel.m_outerRCoordinate = inputParameters.extent[1] / dd4hep::mm;
  paramBarrel.m_outerZCoordinate = inputParameters.extent[3] / dd4hep::mm;
  paramBarrel.m_outerPhiCoordinate = 0.;
  paramBarrel.m_outerSymmetryOrder = 0;
  paramBarrel.m_isMirroredInZ = true;
  paramBarrel.m_nLayers = nlayers;

  paramBarrel.m_layerParametersVector.resize(nlayers);
  float distanceBarrel = inputParameters.extent[0];

  for (unsigned iLayer = 0; iLayer < nlayers; ++iLayer) {
    paramBarrel.m_layerParametersVector.at(iLayer).m_closestDistanceToIp = distanceBarrel / dd4hep::mm;
    paramBarrel.m_layerParametersVector.at(iLayer).m_nRadiationLengths = 0.;
    paramBarrel.m_layerParametersVector.at(iLayer).m_nInteractionLengths = 0.;

    if (layerVecSize > 0)
      distanceBarrel += inputParameters.layers.at(iLayer).sensitive_thickness;
  }

  return;
}

void GeometryCreatorIdea::SetHcalEndcapParameters(const dd4hep::rec::LayeredCalorimeterData& inputParameters,
                                                  PandoraApi::Geometry::SubDetector::Parameters& paramEndcap) const {

  unsigned layerVecSize = 0;
  unsigned nlayers = layerVecSize > 0 ? layerVecSize : 1; // avoid zero

  paramEndcap.m_subDetectorName = "HCalEndcap";
  paramEndcap.m_subDetectorType = pandora::HCAL_ENDCAP;
  paramEndcap.m_innerRCoordinate = inputParameters.extent[0] / dd4hep::mm;
  paramEndcap.m_innerZCoordinate = inputParameters.extent[2] / dd4hep::mm;
  paramEndcap.m_innerPhiCoordinate = 0.;
  paramEndcap.m_innerSymmetryOrder = 0;
  paramEndcap.m_outerRCoordinate = inputParameters.extent[1] / dd4hep::mm;
  paramEndcap.m_outerZCoordinate = inputParameters.extent[3] / dd4hep::mm;
  paramEndcap.m_outerPhiCoordinate = 0.;
  paramEndcap.m_outerSymmetryOrder = 0;
  paramEndcap.m_isMirroredInZ = true;
  paramEndcap.m_nLayers = nlayers;

  paramEndcap.m_layerParametersVector.resize(nlayers);
  float distanceEndcap = inputParameters.extent[2];

  for (unsigned iLayer = 0; iLayer < nlayers; ++iLayer) {
    paramEndcap.m_layerParametersVector.at(iLayer).m_closestDistanceToIp = distanceEndcap / dd4hep::mm;
    paramEndcap.m_layerParametersVector.at(iLayer).m_nRadiationLengths = 0.;
    paramEndcap.m_layerParametersVector.at(iLayer).m_nInteractionLengths = 0.;

    if (layerVecSize > 0)
      distanceEndcap += inputParameters.layers.at(iLayer).sensitive_thickness;
  }

  return;
}
