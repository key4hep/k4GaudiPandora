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
#include "DualReadoutCaloHitCreator.h"

#include "Pandora/PandoraEnumeratedTypes.h"
#include "Pandora/PandoraInputTypes.h"

#include "DDSegmentation/BitFieldCoder.h"

#include "GaudiKernel/AnyDataWrapper.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "k4FWCore/MetadataUtils.h"

DualReadoutCaloHitCreator::DualReadoutCaloHitCreator(const Settings& settings, pandora::Pandora& pandora,
                                                     const Gaudi::Algorithm* algorithm)
    : m_settings(settings), m_pandora(pandora), m_algorithm(*algorithm) {}

pandora::StatusCode DualReadoutCaloHitCreator::createCaloHits(
    const std::vector<std::vector<edm4hep::CalorimeterHit>>& caloHitVectors) const {
  // system decoder to find out which subdetector the hit belongs to
  const dd4hep::DDSegmentation::BitFieldCoder decoderSystem("system:5");

  // loop over hit vectors (one per collection)
  for (size_t iColl = 0; iColl < caloHitVectors.size(); ++iColl) {
    const auto& caloHits = caloHitVectors.at(iColl);

    // no hits in this vector
    if (caloHits.empty())
      continue;

    // loop over calo hits
    for (const auto& hit : caloHits) {
      // check which subdetector the hit belongs to
      const auto& cellID = hit.getCellID();
      uint64_t systemID = decoderSystem.get(cellID, "system");
      // find the corresponding subdetector settings
      auto it = std::find_if(m_settings.m_subDetectorSettings.begin(), m_settings.m_subDetectorSettings.end(),
                             [systemID](const DualReadoutCaloHitCreator::Settings::SubDetectorSettings& settings) {
                               return settings.m_systemID == systemID;
                             });
      if (it == m_settings.m_subDetectorSettings.end()) {
        m_algorithm.warning() << "No subdetector settings found for system ID " << systemID;
        continue;
      }

      const auto& subdetectorSettings = *it;
      const auto& layerFieldName = subdetectorSettings.m_layerFieldName;
      const auto& collectionType = subdetectorSettings.m_collectionType;

      // get collection metadata cellID encoding string
      // const std::string encodingStr =
      //     k4FWCore::getParameter<std::string>(
      //         podio::collMetadataParamName(colName, edm4hep::labels::CellIDEncoding))
      //         .value_or("");
      // const dd4hep::DDSegmentation::BitFieldCoder encoder(encodingStr);
      const std::string encodingStr = subdetectorSettings.m_encodingString;
      const dd4hep::DDSegmentation::BitFieldCoder encoder(encodingStr);

      bool isCherenkov = encoder.get(cellID, m_settings.m_cherenkovFieldName) != 0;
      bool isEcal = collectionType == "ECAL";
      unsigned iLayer = layerFieldName.empty() ? 0 : encoder.get(cellID, layerFieldName);

      // get hit position
      const auto& pos = hit.getPosition();
      float theta = std::atan2(std::sqrt(pos.x * pos.x + pos.y * pos.y), pos.z);
      bool isBarrel = theta > m_settings.m_theta || theta < (M_PI - m_settings.m_theta);

      // create Pandora calo hit
      PandoraApi::CaloHit::Parameters caloHitParameters;
      // see PandoraSDK/include/Pandora/ObjectCreation.h for the full list
      caloHitParameters.m_positionVector = pandora::InputCartesianVector(pandora::CartesianVector(pos.x, pos.y, pos.z));
      caloHitParameters.m_cellGeometry = pandora::InputCellGeometry(pandora::CellGeometry::RECTANGULAR);
      caloHitParameters.m_cellSize0 = subdetectorSettings.m_cellSize; // in mm (only used for PandoraMonitoring)
      caloHitParameters.m_cellSize1 = subdetectorSettings.m_cellSize; // in mm (only used for PandoraMonitoring)
      caloHitParameters.m_time = pandora::InputFloat(hit.getTime());
      caloHitParameters.m_inputEnergy = pandora::InputFloat(hit.getEnergy());
      caloHitParameters.m_hitType = isCherenkov ? pandora::HitType::DRC_CHEREN : pandora::HitType::DRC_SCINT;
      caloHitParameters.m_hitRegion = isBarrel ? pandora::HitRegion::BARREL : pandora::HitRegion::ENDCAP;

      // dummy parameters
      caloHitParameters.m_isDigital = pandora::InputBool(false);
      caloHitParameters.m_layer = layerFieldName.empty() ? 0 : encoder.get(cellID, layerFieldName);
      caloHitParameters.m_expectedDirection = caloHitParameters.m_positionVector.Get().GetUnitVector(); // projective
      caloHitParameters.m_cellNormalVector = caloHitParameters.m_positionVector.Get().GetUnitVector();
      caloHitParameters.m_cellThickness = subdetectorSettings.m_layerThicknesses.empty()
                                              ? subdetectorSettings.m_cellSize
                                              : subdetectorSettings.m_layerThicknesses.at(iLayer);
      caloHitParameters.m_nCellRadiationLengths = 0.;   // not used
      caloHitParameters.m_nCellInteractionLengths = 0.; // not used
      caloHitParameters.m_mipEquivalentEnergy = 0.;     // not used
      // hit type is occupied by the DRC flag, so we use these to distinguish ECAL and HCAL
      caloHitParameters.m_electromagneticEnergy = isEcal ? caloHitParameters.m_inputEnergy : 0.;
      caloHitParameters.m_hadronicEnergy = isEcal ? 0. : caloHitParameters.m_inputEnergy;
      caloHitParameters.m_isInOuterSamplingLayer = false;

      // address to edm4hep calo hit - store address to the hit in the vector
      caloHitParameters.m_pParentAddress = pandora::InputAddress(&hit);

      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
    } // hits
  } // collections

  return pandora::STATUS_CODE_SUCCESS;
}
