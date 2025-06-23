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

#include "Api/PandoraApi.h"

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DDCaloHitCreator.h>
#include <DDRec/DetectorData.h>
#include <DD4hep/DetType.h>
#include <DD4hep/DetectorSelector.h>

#include "GaudiKernel/MsgStream.h"

#include <edm4hep/CalorimeterHit.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// dd4hep::rec::LayeredCalorimeterData * getExtension(std::string detectorName);
dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag = 0);

// double getCoilOuterR();

/// FIXME: HANDLE PROBLEM WHEN EXTENSION IS MISSING
DDCaloHitCreator::DDCaloHitCreator(const Settings& settings, pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm)
    : m_settings(settings), m_pandora(pandora), m_hCalBarrelLayerThickness(0.f), m_hCalEndCapLayerThickness(0.f),
      m_calorimeterHitVector(0), m_volumeManager(), m_algorithm(*algorithm) {
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& barrelLayers =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL),
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))
          ->layers;

  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP),
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))
          ->layers;

  /// Take thicknesses from last layer (was like that before with gear)
  m_hCalEndCapLayerThickness = (endcapLayers.back().inner_thickness + endcapLayers.back().outer_thickness) / dd4hep::mm;
  m_hCalBarrelLayerThickness = (barrelLayers.back().inner_thickness + barrelLayers.back().outer_thickness) / dd4hep::mm;

  if ((m_hCalEndCapLayerThickness < std::numeric_limits<float>::epsilon()) ||
      (m_hCalBarrelLayerThickness < std::numeric_limits<float>::epsilon()))
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
  m_volumeManager = theDetector.volumeManager();
  if (!m_volumeManager.isValid()) {
    theDetector.apply("DD4hepVolumeManager", 0, nullptr);
    m_volumeManager = theDetector.volumeManager();
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDCaloHitCreator::~DDCaloHitCreator() {}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode
DDCaloHitCreator::createCaloHits(const std::map<std::string, std::vector<edm4hep::CalorimeterHit>>& eCaloHitsMap,
                                 const std::vector<edm4hep::CalorimeterHit>& hCalCaloHits,
                                 const std::vector<edm4hep::CalorimeterHit>& muonCaloHits,
                                 const std::vector<edm4hep::CalorimeterHit>& lCalCaloHits,
                                 const std::vector<edm4hep::CalorimeterHit>& lhCalCaloHits) const {
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->createECalCaloHits(eCaloHitsMap))
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->createHCalCaloHits(hCalCaloHits))
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->createMuonCaloHits(muonCaloHits))
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->createLCalCaloHits(lCalCaloHits))
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->createLHCalCaloHits(lhCalCaloHits))

  return pandora::STATUS_CODE_SUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------------------------
 
pandora::StatusCode DDCaloHitCreator::createECalCaloHits(
    const std::map<std::string, std::vector<edm4hep::CalorimeterHit>>& inputECalCaloHits) const {

  std::string initString = "system:5,side:2,module:8,stave:4,layer:9,submodule:4,x:32:-16,y:-16";
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);

  for (const auto& [originalCollectionName, hits] : inputECalCaloHits) {
    if (hits.empty()) {
      continue;
    }
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& barrelLayers =
        getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))
            ->layers;
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers =
        getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))
            ->layers;

    if (barrelLayers.empty() || endcapLayers.empty()) {
      m_algorithm.error() << "Layer information missing for ECAL!" << endmsg;
      return pandora::STATUS_CODE_FAILURE;
    }

    for (const auto& hit : hits) {
      try {

        float absorberCorrection = 1.0;
        float eCalToMip = m_settings.m_eCalToMip;
        float eCalMipThreshold = m_settings.m_eCalMipThreshold;
        float eCalToEMGeV = m_settings.m_eCalToEMGeV;
        float eCalToHadGeVBarrel = m_settings.m_eCalToHadGeVBarrel;
        float eCalToHadGeVEndCap = m_settings.m_eCalToHadGeVEndCap;

        if (m_settings.m_useEcalScLayers) {
          std::string collectionName = originalCollectionName;
          std::ranges::transform(collectionName, collectionName.begin(),
                                 [](unsigned char c) { return std::tolower(c); });

          if (originalCollectionName.find("ecal") == std::string::npos) {
            m_algorithm.warning() << "WARNING: Mismatching hybrid ECal collection name: " << originalCollectionName << endmsg;
          }

          if (originalCollectionName.find("si") != std::string::npos) {
            eCalToMip = m_settings.m_eCalSiToMip;
            eCalMipThreshold = m_settings.m_eCalSiMipThreshold;
            eCalToEMGeV = m_settings.m_eCalSiToEMGeV;
            eCalToHadGeVBarrel = m_settings.m_eCalSiToHadGeVBarrel;
            eCalToHadGeVEndCap = m_settings.m_eCalSiToHadGeVEndCap;
          } else if (originalCollectionName.find("sc") != std::string::npos) {
            eCalToMip = m_settings.m_eCalScToMip;
            eCalMipThreshold = m_settings.m_eCalScMipThreshold;
            eCalToEMGeV = m_settings.m_eCalScToEMGeV;
            eCalToHadGeVBarrel = m_settings.m_eCalScToHadGeVBarrel;
            eCalToHadGeVEndCap = m_settings.m_eCalScToHadGeVEndCap;
          }
        }

        PandoraApi::CaloHit::Parameters caloHitParameters;
        caloHitParameters.m_hitType = pandora::ECAL;
        caloHitParameters.m_isDigital = false;
        caloHitParameters.m_layer = bitFieldCoder.get(hit.getCellID(), "layer");
        caloHitParameters.m_isInOuterSamplingLayer = false;
        this->getCommonCaloHitProperties(hit, caloHitParameters);

        if ((!m_settings.m_useSystemId && std::fabs(hit.getPosition()[2]) < m_settings.m_eCalBarrelOuterZ) ||
            (m_settings.m_useSystemId && bitFieldCoder.get(hit.getCellID(), "system") == m_settings.m_ecalBarrelSystemId)) {
          this->getBarrelCaloHitProperties(hit, barrelLayers, m_settings.m_eCalBarrelInnerSymmetry, caloHitParameters,
                                           m_settings.m_eCalBarrelNormalVector, absorberCorrection);
          caloHitParameters.m_hadronicEnergy = eCalToHadGeVBarrel * hit.getEnergy();
        } else {
          this->getEndCapCaloHitProperties(hit, endcapLayers, caloHitParameters, absorberCorrection);
          caloHitParameters.m_hadronicEnergy = eCalToHadGeVEndCap * hit.getEnergy();
        }

        caloHitParameters.m_mipEquivalentEnergy = hit.getEnergy() * eCalToMip * absorberCorrection;

        if (caloHitParameters.m_mipEquivalentEnergy.Get() < eCalMipThreshold)
          continue;

        caloHitParameters.m_electromagneticEnergy = eCalToEMGeV * hit.getEnergy();

        if (m_settings.m_stripSplittingOn) {
          const float splitCellSize =
              std::min(caloHitParameters.m_cellSize0.Get(), caloHitParameters.m_cellSize1.Get());
          caloHitParameters.m_cellSize0 = splitCellSize;
          caloHitParameters.m_cellSize1 = splitCellSize;
        }

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                PandoraApi::CaloHit::Create(m_pandora, caloHitParameters))

      } catch (const std::exception& e) {
        m_algorithm.error() << "Exception processing ECAL hit: " << e.what() << endmsg;
      }
    }
  }
  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreator::createHCalCaloHits(const std::vector<edm4hep::CalorimeterHit>& hCalCaloHits) const {

  if (hCalCaloHits.empty())
    return pandora::STATUS_CODE_SUCCESS;

  // Retrieve calorimeter layers from dd4hep
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& barrelLayers =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::BARREL),
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))
          ->layers;
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::HADRONIC | dd4hep::DetType::ENDCAP),
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))
          ->layers;

  if (barrelLayers.empty() || endcapLayers.empty()) {
    m_algorithm.error() << "Layer information missing for HCAL!" << endmsg;
    return pandora::STATUS_CODE_FAILURE;
  }

  std::string initString = "system:5,side:2,module:8,stave:4,layer:9,submodule:4,x:32:-16,y:-16";
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);

  for (const auto& hit : hCalCaloHits) {
    try {
      PandoraApi::CaloHit::Parameters caloHitParameters;
      caloHitParameters.m_hitType = pandora::HCAL;
      caloHitParameters.m_isDigital = false;
      caloHitParameters.m_layer = bitFieldCoder.get(hit.getCellID(), "layer");
      caloHitParameters.m_isInOuterSamplingLayer = (this->getNLayersFromEdge(hit) <= m_settings.m_nOuterSamplingLayers);
      this->getCommonCaloHitProperties(hit, caloHitParameters);

      float absorberCorrection = 1.0;

      if (std::fabs(hit.getPosition()[2]) < m_settings.m_hCalBarrelOuterZ) {
        this->getBarrelCaloHitProperties(hit, barrelLayers, m_settings.m_hCalBarrelInnerSymmetry, caloHitParameters,
                                         m_settings.m_hCalBarrelNormalVector, absorberCorrection);
      } else {
        this->getEndCapCaloHitProperties(hit, endcapLayers, caloHitParameters, absorberCorrection);
      }

      caloHitParameters.m_mipEquivalentEnergy = hit.getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

      if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
        continue;

      caloHitParameters.m_hadronicEnergy =
          std::min(m_settings.m_hCalToHadGeV * hit.getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
      caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * hit.getEnergy();

      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::CaloHit::Create(m_pandora, caloHitParameters))

    } catch (const std::exception& e) {
      m_algorithm.error() << "Exception processing HCAL hit: " << e.what() << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreator::createMuonCaloHits(const std::vector<edm4hep::CalorimeterHit>& muonCaloHits) const {

  if (muonCaloHits.empty())
    return pandora::STATUS_CODE_SUCCESS;

  // Initialize BitFieldCoder with proper encoding string
  std::string initString = "system:5,side:2,module:8,stave:4,layer:9,submodule:4,x:32:-16,y:-16";
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);

  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& barrelLayers =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::BARREL),
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))
          ->layers;
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers =
      getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::MUON | dd4hep::DetType::ENDCAP),
                   (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD))
          ->layers;

  if (barrelLayers.empty() || endcapLayers.empty()) {
    m_algorithm.error() << "Layer information missing for MUON!" << endmsg;
    return pandora::STATUS_CODE_FAILURE;
  }

  for (const auto& hit : muonCaloHits) {
    try {
      PandoraApi::CaloHit::Parameters caloHitParameters;
      caloHitParameters.m_hitType = pandora::MUON;
      caloHitParameters.m_layer = bitFieldCoder.get(hit.getCellID(), "layer");
      caloHitParameters.m_isInOuterSamplingLayer = true;
      this->getCommonCaloHitProperties(hit, caloHitParameters);

      const float radius =
          std::sqrt(hit.getPosition()[0] * hit.getPosition()[0] + hit.getPosition()[1] * hit.getPosition()[1]);

      const bool isWithinCoil = (radius < m_settings.m_coilOuterR);
      const bool isInBarrelRegion = (std::fabs(hit.getPosition()[2]) < m_settings.m_muonBarrelOuterZ);

      float absorberCorrection = 1.0;

      if (isInBarrelRegion && isWithinCoil) {
        m_algorithm.warning() << "BIG WARNING: CANNOT HANDLE PLUG HITS (no plug in CLIC model), DO NOTHING!" << endmsg;
      } else if (isInBarrelRegion) {
        this->getBarrelCaloHitProperties(hit, barrelLayers, m_settings.m_muonBarrelInnerSymmetry, caloHitParameters,
                                         m_settings.m_muonBarrelNormalVector, absorberCorrection);
      } else {
        this->getEndCapCaloHitProperties(hit, endcapLayers, caloHitParameters, absorberCorrection);
      }

      if (m_settings.m_muonDigitalHits > 0) {
        caloHitParameters.m_isDigital = true;
        caloHitParameters.m_inputEnergy = m_settings.m_muonHitEnergy;
        caloHitParameters.m_hadronicEnergy = m_settings.m_muonHitEnergy;
        caloHitParameters.m_electromagneticEnergy = m_settings.m_muonHitEnergy;
        caloHitParameters.m_mipEquivalentEnergy = 1.0f;
      } else {
        caloHitParameters.m_isDigital = false;
        caloHitParameters.m_inputEnergy = hit.getEnergy();
        caloHitParameters.m_hadronicEnergy = hit.getEnergy();
        caloHitParameters.m_electromagneticEnergy = hit.getEnergy();
        caloHitParameters.m_mipEquivalentEnergy = hit.getEnergy() * m_settings.m_muonToMip;
      }

      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::CaloHit::Create(m_pandora, caloHitParameters))

    } catch (const std::exception& e) {
      m_algorithm.error() << "Exception processing MUON hit: " << e.what() << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode
DDCaloHitCreator::createLCalCaloHits(const std::vector<edm4hep::CalorimeterHit>& inputLCalCaloHits) const {

  if (inputLCalCaloHits.empty())
    return pandora::STATUS_CODE_SUCCESS;

  std::string initString = "system:5,side:2,module:8,stave:4,layer:9,submodule:4,x:32:-16,y:-16";
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);

  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers =
      getExtension(dd4hep::DetType::CALORIMETER | dd4hep::DetType::ENDCAP | dd4hep::DetType::ELECTROMAGNETIC |
                       dd4hep::DetType::FORWARD,
                   dd4hep::DetType::AUXILIARY)
          ->layers;

  if (endcapLayers.empty()) {
    m_algorithm.error() << "Layer information missing for LCal!" << endmsg;
    return pandora::STATUS_CODE_FAILURE;
  }

  for (const auto& hit : inputLCalCaloHits) {
    try {
      PandoraApi::CaloHit::Parameters caloHitParameters;
      caloHitParameters.m_hitType = pandora::ECAL;
      caloHitParameters.m_isDigital = false;
      caloHitParameters.m_layer = bitFieldCoder.get(hit.getCellID(), "layer");
      caloHitParameters.m_isInOuterSamplingLayer = false;
      this->getCommonCaloHitProperties(hit, caloHitParameters);

      float absorberCorrection = 1.0;
      this->getEndCapCaloHitProperties(hit, endcapLayers, caloHitParameters, absorberCorrection);

      caloHitParameters.m_mipEquivalentEnergy = hit.getEnergy() * m_settings.m_eCalToMip * absorberCorrection;

      if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
        continue;

      caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * hit.getEnergy();
      caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVEndCap * hit.getEnergy();

      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::CaloHit::Create(m_pandora, caloHitParameters))

    } catch (const std::exception& e) {
      m_algorithm.error() << "Exception processing LCal hit: " << e.what() << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDCaloHitCreator::createLHCalCaloHits(const std::vector<edm4hep::CalorimeterHit>& LHCalCaloHits) const {

  if (LHCalCaloHits.empty())
    return pandora::STATUS_CODE_SUCCESS;

  std::string initString = "system:5,side:2,module:8,stave:4,layer:9,submodule:4,x:32:-16,y:-16";
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);

  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& endcapLayers =
      getExtension(dd4hep::DetType::CALORIMETER | dd4hep::DetType::ENDCAP | dd4hep::DetType::HADRONIC |
                       dd4hep::DetType::FORWARD,
                   dd4hep::DetType::AUXILIARY)
          ->layers;

  if (endcapLayers.empty()) {
    m_algorithm.error() << "Layer information missing for LHCal!" << endmsg;
    return pandora::STATUS_CODE_FAILURE;
  }

  for (const auto& hit : LHCalCaloHits) {
    try {
      PandoraApi::CaloHit::Parameters caloHitParameters;
      caloHitParameters.m_hitType = pandora::HCAL;
      caloHitParameters.m_isDigital = false;
      caloHitParameters.m_layer = bitFieldCoder.get(hit.getCellID(), "layer");
      caloHitParameters.m_isInOuterSamplingLayer = (this->getNLayersFromEdge(hit) <= m_settings.m_nOuterSamplingLayers);
      this->getCommonCaloHitProperties(hit, caloHitParameters);

      float absorberCorrection = 1.0;
      this->getEndCapCaloHitProperties(hit, endcapLayers, caloHitParameters, absorberCorrection);

      caloHitParameters.m_mipEquivalentEnergy = hit.getEnergy() * m_settings.m_hCalToMip * absorberCorrection;

      if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_hCalMipThreshold)
        continue;

      caloHitParameters.m_hadronicEnergy =
          std::min(m_settings.m_hCalToHadGeV * hit.getEnergy(), m_settings.m_maxHCalHitHadronicEnergy);
      caloHitParameters.m_electromagneticEnergy = m_settings.m_hCalToEMGeV * hit.getEnergy();

      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::CaloHit::Create(m_pandora, caloHitParameters))

    } catch (const std::exception& e) {
      m_algorithm.error() << "Exception processing LHCal hit: " << e.what() << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreator::getCommonCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                                  PandoraApi::CaloHit::Parameters& caloHitParameters) const {
  const auto position = hit.getPosition();
  const pandora::CartesianVector positionVector(position.x, position.y, position.z);

  caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
  caloHitParameters.m_positionVector = positionVector;
  caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
  caloHitParameters.m_pParentAddress = static_cast<const void*>(&hit);
  caloHitParameters.m_inputEnergy = hit.getEnergy();
  caloHitParameters.m_time = hit.getTime();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreator::getEndCapCaloHitProperties(
    const edm4hep::CalorimeterHit& hit, const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
    PandoraApi::CaloHit::Parameters& caloHitParameters, float& absorberCorrection) const {
  caloHitParameters.m_hitRegion = pandora::ENDCAP;

  // Determine physical layer index
  const int physicalLayer =
      std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size() - 1));

  caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0 / dd4hep::mm;
  caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1 / dd4hep::mm;

  double thickness =
      (layers[physicalLayer].inner_thickness + layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
  double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
  double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;
  double layerAbsorberThickness =
      (layers[physicalLayer].inner_thickness - layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;

  if (physicalLayer > 0) {
    thickness +=
        (layers[physicalLayer - 1].outer_thickness - layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
    nRadLengths += layers[physicalLayer - 1].outer_nRadiationLengths;
    nIntLengths += layers[physicalLayer - 1].outer_nInteractionLengths;
    layerAbsorberThickness +=
        (layers[physicalLayer - 1].outer_thickness - layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
  }

  caloHitParameters.m_cellThickness = thickness;
  caloHitParameters.m_nCellRadiationLengths = nRadLengths;
  caloHitParameters.m_nCellInteractionLengths = nIntLengths;

  if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() ||
      caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon()) {
    m_algorithm.warning() << "CaloHitCreator::GetEndCapCaloHitProperties: Calo hit has 0 radiation length or interaction length. Not creating a Pandora calo hit." << endmsg;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }

  // Compute absorber correction factor
  absorberCorrection = 1.0;
  for (size_t i = 0; i < layers.size(); ++i) {
    float absorberThickness = (layers[i].inner_thickness - layers[i].sensitive_thickness / 2.0) / dd4hep::mm;

    if (i > 0) {
      absorberThickness += (layers[i - 1].outer_thickness - layers[i - 1].sensitive_thickness / 2.0) / dd4hep::mm;
    }

    if (absorberThickness < std::numeric_limits<float>::epsilon()) {
      continue;
    }

    if (layerAbsorberThickness > std::numeric_limits<float>::epsilon()) {
      absorberCorrection = absorberThickness / layerAbsorberThickness;
    }

    break;
  }

  // Set normal vector for endcap hits
  caloHitParameters.m_cellNormalVector =
      (hit.getPosition()[2] > 0) ? pandora::CartesianVector(0, 0, 1) : pandora::CartesianVector(0, 0, -1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDCaloHitCreator::getBarrelCaloHitProperties(
    const edm4hep::CalorimeterHit& hit, const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
    unsigned int barrelSymmetryOrder, PandoraApi::CaloHit::Parameters& caloHitParameters,
    const std::vector<float>& normalVector, float& absorberCorrection) const {
  caloHitParameters.m_hitRegion = pandora::BARREL;

  const int physicalLayer =
      std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size() - 1));
  caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0 / dd4hep::mm;
  caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1 / dd4hep::mm;

  double thickness =
      (layers[physicalLayer].inner_thickness + layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
  double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
  double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;
  double layerAbsorberThickness =
      (layers[physicalLayer].inner_thickness - layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;

  if (physicalLayer > 0) {
    thickness +=
        (layers[physicalLayer - 1].outer_thickness - layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
    nRadLengths += layers[physicalLayer - 1].outer_nRadiationLengths;
    nIntLengths += layers[physicalLayer - 1].outer_nInteractionLengths;
    layerAbsorberThickness +=
        (layers[physicalLayer - 1].outer_thickness - layers[physicalLayer].sensitive_thickness / 2.0) / dd4hep::mm;
  }

  caloHitParameters.m_cellThickness = thickness;
  caloHitParameters.m_nCellRadiationLengths = nRadLengths;
  caloHitParameters.m_nCellInteractionLengths = nIntLengths;

  if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() ||
      caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon()) {
    m_algorithm.warning() << "CaloHitCreator::getBarrelCaloHitProperties: Calo hit has 0 radiation or interaction length." << endmsg;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }

  absorberCorrection = 1.0;
  for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i) {
    float absorberThickness = (layers[i].inner_thickness - layers[i].sensitive_thickness / 2.0) / dd4hep::mm;
    if (i > 0)
      absorberThickness += (layers[i - 1].outer_thickness - layers[i - 1].sensitive_thickness / 2.0) / dd4hep::mm;
    if (absorberThickness < std::numeric_limits<float>::epsilon())
      continue;
    if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
      absorberCorrection = absorberThickness / layerAbsorberThickness;
    break;
  }

  if (barrelSymmetryOrder > 2) {
    if (hit.getCellID() != 0) {
      auto staveDetElement = m_volumeManager.lookupDetElement(hit.getCellID());
      dd4hep::Position local1(0.0, 0.0, 0.0);
      dd4hep::Position local2(normalVector[0], normalVector[1], normalVector[2]);
      dd4hep::Position global1(0.0, 0.0, 0.0);
      dd4hep::Position global2(0.0, 0.0, 0.0);
      staveDetElement.nominal().localToWorld(local1, global1);
      staveDetElement.nominal().localToWorld(local2, global2);
      dd4hep::Position normal(global2 - global1);

      m_algorithm.debug() << "DetElement: " << staveDetElement.name()
                          << " Parent: " << staveDetElement.parent().name()
                          << " Grandparent: " << staveDetElement.parent().parent().name()
                          << " CellID: " << hit.getCellID()
                          << " PhiLoc: " << atan2(global1.y(), global1.x()) * 180 / M_PI
                          << " PhiNor: " << atan2(normal.y(), normal.x()) * 180 / M_PI
                          << " Normal Vector: " << normal.x() << " " << normal.y() << " " << normal.z()
                          << endmsg;

      caloHitParameters.m_cellNormalVector = pandora::CartesianVector(normal.x(), normal.y(), normal.z());
    } else {
      const double phi = atan2(hit.getPosition()[1], hit.getPosition()[0]);
      m_algorithm.warning() << "Hit does not have any cellIDs set, using phi-direction for normal vector "
                            << " Phi: " << phi * 180 / M_PI << endmsg;
      caloHitParameters.m_cellNormalVector = pandora::CartesianVector(std::cos(phi), std::sin(phi), 0.0);
    }
  } else {
    const double phi = atan2(hit.getPosition()[1], hit.getPosition()[0]);
    caloHitParameters.m_cellNormalVector = pandora::CartesianVector(std::cos(phi), std::sin(phi), 0);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DDCaloHitCreator::getNLayersFromEdge(const edm4hep::CalorimeterHit& hit) const {
  // Calo hit coordinate calculations
  const float barrelMaximumRadius =
      this->getMaximumRadius(hit, m_settings.m_hCalBarrelOuterSymmetry, m_settings.m_hCalBarrelOuterPhi0);
  const float endCapMaximumRadius =
      this->getMaximumRadius(hit, m_settings.m_hCalEndCapInnerSymmetryOrder, m_settings.m_hCalEndCapInnerPhiCoordinate);
  const float caloHitAbsZ = std::fabs(hit.getPosition()[2]);

  // Distance from radial outer
  float radialDistanceToEdge = std::numeric_limits<float>::max();

  if (caloHitAbsZ < m_settings.m_eCalBarrelOuterZ) {
    radialDistanceToEdge = (m_settings.m_hCalBarrelOuterR - barrelMaximumRadius) / m_hCalBarrelLayerThickness;
  } else {
    radialDistanceToEdge = (m_settings.m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness;
  }

  // Distance from rear of endcap outer
  float rearDistanceToEdge = std::numeric_limits<float>::max();

  if (caloHitAbsZ >= m_settings.m_eCalBarrelOuterZ) {
    rearDistanceToEdge = (m_settings.m_hCalEndCapOuterZ - caloHitAbsZ) / m_hCalEndCapLayerThickness;
  } else {
    const float rearDistance = (m_settings.m_eCalBarrelOuterZ - caloHitAbsZ) / m_hCalBarrelLayerThickness;

    if (rearDistance < m_settings.m_layersFromEdgeMaxRearDistance) {
      const float overlapDistance = (m_settings.m_hCalEndCapOuterR - endCapMaximumRadius) / m_hCalEndCapLayerThickness;
      rearDistanceToEdge = std::max(rearDistance, overlapDistance);
    }
  }

  return static_cast<int>(std::min(radialDistanceToEdge, rearDistanceToEdge));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DDCaloHitCreator::getMaximumRadius(const edm4hep::CalorimeterHit& pCaloHit, const unsigned int symmetryOrder,
                                         const float phi0) const {
  const auto& pCaloHitPosition(pCaloHit.getPosition());

  if (symmetryOrder <= 2)
    return std::sqrt((pCaloHitPosition[0] * pCaloHitPosition[0]) + (pCaloHitPosition[1] * pCaloHitPosition[1]));

  float maximumRadius = 0.f;
  const float twoPi = 2.f * M_PI;

  for (unsigned int i = 0; i < symmetryOrder; ++i) {
    const float phi = phi0 + i * twoPi / static_cast<float>(symmetryOrder);
    float radius = pCaloHitPosition[0] * std::cos(phi) + pCaloHitPosition[1] * std::sin(phi);

    if (radius > maximumRadius)
      maximumRadius = radius;
  }

  return maximumRadius;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDCaloHitCreator::Settings::Settings()
    : m_eCalCaloHitCollections(StringVector()), m_hCalCaloHitCollections(StringVector()),
      m_lCalCaloHitCollections(StringVector()), m_lHCalCaloHitCollections(StringVector()),
      m_muonCaloHitCollections(StringVector()), m_eCalToMip(1.f), m_hCalToMip(1.f), m_muonToMip(1.f),
      m_eCalMipThreshold(0.f), m_hCalMipThreshold(0.f), m_muonMipThreshold(0.f), m_eCalToEMGeV(1.f),
      m_eCalToHadGeVBarrel(1.f), m_eCalToHadGeVEndCap(1.f), m_hCalToEMGeV(1.f), m_hCalToHadGeV(1.f),
      m_muonDigitalHits(1), m_muonHitEnergy(0.5f), m_maxHCalHitHadronicEnergy(10000.f), m_nOuterSamplingLayers(3),
      m_layersFromEdgeMaxRearDistance(250.f), m_hCalEndCapInnerSymmetryOrder(4), m_hCalEndCapInnerPhiCoordinate(0.f),
      m_stripSplittingOn(0), m_useEcalScLayers(0), m_eCalSiToMip(1.f), m_eCalScToMip(1.f), m_eCalSiMipThreshold(0.f),
      m_eCalScMipThreshold(0.f), m_eCalSiToEMGeV(1.f), m_eCalScToEMGeV(1.f), m_eCalSiToHadGeVBarrel(1.f),
      m_eCalScToHadGeVBarrel(1.f), m_eCalSiToHadGeVEndCap(1.f), m_eCalScToHadGeVEndCap(1.f), m_eCalBarrelOuterZ(0.f),
      m_hCalBarrelOuterZ(0.f), m_muonBarrelOuterZ(0.f), m_coilOuterR(0.f), m_eCalBarrelInnerPhi0(0.f),
      m_eCalBarrelInnerSymmetry(0.f), m_hCalBarrelInnerPhi0(0.f), m_hCalBarrelInnerSymmetry(0.f),
      m_muonBarrelInnerPhi0(0.f), m_muonBarrelInnerSymmetry(0.f), m_hCalEndCapOuterR(0.f), m_hCalEndCapOuterZ(0.f),
      m_hCalBarrelOuterR(0.f), m_hCalBarrelOuterPhi0(0.f), m_hCalBarrelOuterSymmetry(0.f), m_useSystemId(false),
      m_ecalBarrelSystemId(-1), m_hcalBarrelSystemId(-1), m_eCalBarrelNormalVector({0.0, 0.0, 1.0}),
      m_hCalBarrelNormalVector({0.0, 0.0, 1.0}), m_muonBarrelNormalVector({0.0, 0.0, 1.0}) {}
