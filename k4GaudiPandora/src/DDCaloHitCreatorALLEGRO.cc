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

#include "DDCaloHitCreatorALLEGRO.h"
#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetElement.h>
#include <DD4hep/DetType.h>
#include <DD4hep/Detector.h>
#include <DD4hep/DetectorSelector.h>
#include <DDRec/DetectorData.h>

DDCaloHitCreatorALLEGRO::DDCaloHitCreatorALLEGRO(const Settings& settings, pandora::Pandora& pandora,
                                                 const Gaudi::Algorithm* algorithm)
    : DDCaloHitCreator(settings, pandora, algorithm), m_hcalBarrelSegmentation(nullptr),
      m_hcalEndcapSegmentation(nullptr) {

  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();

  const std::vector<dd4hep::DetElement>& hCalBarrelDetector =
      dd4hep::DetectorSelector(mainDetector)
          .detectors((dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

  const std::vector<dd4hep::DetElement>& hCalEndcapDetector =
      dd4hep::DetectorSelector(mainDetector)
          .detectors((dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

  auto hCalBarrelSensDet = mainDetector.sensitiveDetector(hCalBarrelDetector.at(0).name());
  if (hCalBarrelSensDet.isValid()) {
    auto readout = hCalBarrelSensDet.readout();
    if (readout.isValid()) {
      m_algorithm.debug() << "ALLEGRO DDCaloHitCreator: Readout: " << readout.name() << endmsg;
      // get segmentation
      dd4hep::DDSegmentation::Segmentation* aSegmentation = readout.segmentation().segmentation();
      if (aSegmentation == nullptr) {
        m_algorithm.error() << "ALLEGRO DDCaloHitCreator: Segmentation does not exist." << endmsg;
      }

      std::string segmentationType = aSegmentation->type();
      m_algorithm.debug() << "ALLEGRO DDCaloHitCreator: Segmentation type : " << segmentationType << endmsg;
      m_hcalBarrelSegmentation = std::shared_ptr<dd4hep::DDSegmentation::Segmentation>(
          aSegmentation, [](dd4hep::DDSegmentation::Segmentation*) {});
    }
  }

  auto hCalEndcapSensDet = mainDetector.sensitiveDetector(hCalEndcapDetector.at(0).name());
  if (hCalEndcapSensDet.isValid()) {
    auto readout = hCalEndcapSensDet.readout();
    if (readout.isValid()) {
      m_algorithm.debug() << "ALLEGRO DDCaloHitCreator: Readout: " << readout.name() << endmsg;
      // get segmentation
      dd4hep::DDSegmentation::Segmentation* aSegmentation = readout.segmentation().segmentation();
      if (aSegmentation == nullptr) {
        m_algorithm.error() << "ALLEGRO DDCaloHitCreator: Segmentation does not exist." << endmsg;
      }

      std::string segmentationType = aSegmentation->type();
      m_algorithm.debug() << "ALLEGRO DDCaloHitCreator: Segmentation type : " << segmentationType << endmsg;
      m_hcalEndcapSegmentation = std::shared_ptr<dd4hep::DDSegmentation::Segmentation>(
          aSegmentation, [](dd4hep::DDSegmentation::Segmentation*) {});
    }
  }

  // make sure that we have got readout segmentation objects
  if (!m_hcalBarrelSegmentation || !m_hcalEndcapSegmentation) {
    m_algorithm.error() << "ALLEGRO DDCaloHitCreator: Unable to get segmentation objects." << endmsg;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }
}

pandora::StatusCode DDCaloHitCreatorALLEGRO::createCaloHits(const std::vector<CollectionDescriptor>& eCalCaloHits,
                                                            const std::vector<CollectionDescriptor>& hCalCaloHits,
                                                            const std::vector<CollectionDescriptor>& muonCaloHits,
                                                            const std::vector<edm4hep::CalorimeterHit>&,
                                                            const std::vector<edm4hep::CalorimeterHit>&) const {
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, createECalCaloHits(eCalCaloHits))
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, createHCalCaloHits(hCalCaloHits))
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, createMuonCaloHits(muonCaloHits))

  return pandora::STATUS_CODE_SUCCESS;
}

void DDCaloHitCreatorALLEGRO::getCommonCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                                         PandoraApi::CaloHit::Parameters& caloHitParameters) const {
  const auto position = hit.getPosition();
  const pandora::CartesianVector positionVector(position.x, position.y, position.z);

  // AD: HCAL cells are rectangular. If the hit is in the HCAL get the dimensions from the segmentation class
  if (caloHitParameters.m_hitType.Get() == pandora::HCAL)
    caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
  else
    caloHitParameters.m_cellGeometry = pandora::POINTING;

  caloHitParameters.m_positionVector = positionVector;
  caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
  caloHitParameters.m_pParentAddress = static_cast<const void*>(&hit);
  caloHitParameters.m_inputEnergy = hit.getEnergy();
  caloHitParameters.m_time = hit.getTime();
}

void DDCaloHitCreatorALLEGRO::getEndCapCaloHitProperties(
    const edm4hep::CalorimeterHit& hit, const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
    PandoraApi::CaloHit::Parameters& caloHitParameters, float&) const {
  caloHitParameters.m_hitRegion = pandora::ENDCAP;

  m_algorithm.debug() << "ALLEGRO getEndCapCaloHitProperties:" << endmsg;
  m_algorithm.debug() << "  Hit type: " << caloHitParameters.m_hitType.Get() << endmsg;
  m_algorithm.debug() << "  Hit layer: " << caloHitParameters.m_layer.Get() << endmsg;
  m_algorithm.debug() << "  Hit cellGeometry: " << caloHitParameters.m_cellGeometry.Get() << endmsg;

  const int physicalLayer(
      std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size() - 1)));

  // AD: HCAL cells are rectangular. if the hit is in the HCAL get the dimensions from the segmentation class
  if (caloHitParameters.m_hitType.Get() == pandora::HCAL) {
    std::vector<double> cellSize = m_hcalEndcapSegmentation->cellDimensions(hit.getCellID());
    caloHitParameters.m_cellSize0 = cellSize[0] / dd4hep::mm;
    caloHitParameters.m_cellSize1 = cellSize[1] / dd4hep::mm;
  } else {
    caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0;
    caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1;
  }

  m_algorithm.debug() << "  Hit cellSize0: " << caloHitParameters.m_cellSize0.Get() << endmsg;
  m_algorithm.debug() << "  Hit cellSize1: " << caloHitParameters.m_cellSize1.Get() << endmsg;
  m_algorithm.debug() << "  Nlayers: " << layers.size() << endmsg;

  double thickness =
      (layers[physicalLayer].inner_thickness + layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
  double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
  double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;

  if (physicalLayer > 0) {
    thickness +=
        (layers[physicalLayer - 1].outer_thickness - layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
    nRadLengths += layers[physicalLayer - 1].outer_nRadiationLengths;
    nIntLengths += layers[physicalLayer - 1].outer_nInteractionLengths;
  }

  caloHitParameters.m_cellThickness = thickness;
  caloHitParameters.m_nCellRadiationLengths = nRadLengths;
  caloHitParameters.m_nCellInteractionLengths = nIntLengths;

  if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() ||
      caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon()) {
    m_algorithm.warning()
        << "ALLEGRO CaloHitCreator::getEndCapCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit."
        << endmsg;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }
  caloHitParameters.m_cellNormalVector =
      (hit.getPosition()[2] > 0) ? pandora::CartesianVector(0, 0, 1) : pandora::CartesianVector(0, 0, -1);

  m_algorithm.debug() << "ALLEGRO getEndCapCaloHitProperties: physLayer: " << physicalLayer
                      << " layer: " << caloHitParameters.m_layer.Get()
                      << " nX0: " << caloHitParameters.m_nCellRadiationLengths.Get()
                      << " nLambdaI: " << caloHitParameters.m_nCellInteractionLengths.Get()
                      << " thickness: " << caloHitParameters.m_cellThickness.Get() << endmsg;
}

//------------------------------------------------------------------------------------------------------------------------------------------
void DDCaloHitCreatorALLEGRO::getBarrelCaloHitProperties(
    const edm4hep::CalorimeterHit& hit, const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
    unsigned int, PandoraApi::CaloHit::Parameters& caloHitParameters, FloatVector const&, float&) const {
  caloHitParameters.m_hitRegion = pandora::BARREL;

  m_algorithm.debug() << "ALLEGRO getBarrelCaloHitProperties:" << endmsg;
  m_algorithm.debug() << "  Hit type: " << caloHitParameters.m_hitType.Get() << endmsg;
  m_algorithm.debug() << "  Hit layer: " << caloHitParameters.m_layer.Get() << endmsg;
  m_algorithm.debug() << "  Hit cellGeometry: " << caloHitParameters.m_cellGeometry.Get() << endmsg;

  const int physicalLayer(
      std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size() - 1)));

  // AD: HCAL cells are rectangular. If the hit is in the HCAL get dimensions from the segmentation class
  if (caloHitParameters.m_hitType.Get() == pandora::HCAL) {
    std::vector<double> cellSize = m_hcalBarrelSegmentation->cellDimensions(hit.getCellID());
    caloHitParameters.m_cellSize0 = cellSize[0] / dd4hep::mm;
    caloHitParameters.m_cellSize1 = cellSize[1] / dd4hep::mm;
  } else {
    caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0;
    caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1;
  }

  m_algorithm.debug() << "  Hit cellSize0: " << caloHitParameters.m_cellSize0.Get() << endmsg;
  m_algorithm.debug() << "  Hit cellSize1: " << caloHitParameters.m_cellSize1.Get() << endmsg;
  m_algorithm.debug() << "  Nlayers: " << layers.size() << endmsg;

  double thickness =
      (layers[physicalLayer].inner_thickness + layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
  double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
  double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;

  if (physicalLayer > 0) {
    thickness +=
        (layers[physicalLayer - 1].outer_thickness - layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
    nRadLengths += layers[physicalLayer - 1].outer_nRadiationLengths;
    nIntLengths += layers[physicalLayer - 1].outer_nInteractionLengths;
  }

  caloHitParameters.m_cellThickness = thickness;
  caloHitParameters.m_nCellRadiationLengths = nRadLengths;
  caloHitParameters.m_nCellInteractionLengths = nIntLengths;

  if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() ||
      caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon()) {
    m_algorithm.warning()
        << "ALLEGRO CaloHitCreator::getBarrelCaloHitProperties Calo hit has 0 radiation length or interaction length: \
            not creating a Pandora calo hit."
        << endmsg;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }

  const auto position(hit.getPosition());
  const float phi = std::atan2(position.y, position.y);
  caloHitParameters.m_cellNormalVector = pandora::CartesianVector(std::cos(phi), std::sin(phi), 0);

  m_algorithm.debug() << "ALLEGRO getBarrelCaloHitProperties: physLayer: " << physicalLayer
                      << " layer: " << caloHitParameters.m_layer.Get()
                      << " nX0: " << caloHitParameters.m_nCellRadiationLengths.Get()
                      << " nLambdaI: " << caloHitParameters.m_nCellInteractionLengths.Get()
                      << " thickness: " << caloHitParameters.m_cellThickness.Get() << endmsg;
}
