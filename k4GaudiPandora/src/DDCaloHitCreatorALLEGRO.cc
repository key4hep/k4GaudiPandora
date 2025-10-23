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
 *  @file   DDMarlinPandora/src/DDCaloHitCreator.cc
 *
 *  @brief  Implementation of the calo hit creator class.
 *
 *  $Log: $
 */

#include "DDCaloHitCreatorALLEGRO.h"

#include <DDRec/DetectorData.h>

// dd4hep::rec::LayeredCalorimeterData * getExtension(std::string detectorName);
dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag = 0);

DDCaloHitCreatorALLEGRO::DDCaloHitCreatorALLEGRO(const Settings& settings, pandora::Pandora& pandora,
                                                 const Gaudi::Algorithm* algorithm)
    : DDCaloHitCreator(settings, pandora, algorithm) {}

pandora::StatusCode DDCaloHitCreatorALLEGRO::createCaloHits(
    const std::map<std::string, std::vector<edm4hep::CalorimeterHit>>& eCaloHitsMap,
    const std::vector<edm4hep::CalorimeterHit>& hCalCaloHits, const std::vector<edm4hep::CalorimeterHit>& muonCaloHits,
    const std::vector<edm4hep::CalorimeterHit>&, const std::vector<edm4hep::CalorimeterHit>&) const {
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, createECalCaloHits(eCaloHitsMap))
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, createHCalCaloHits(hCalCaloHits))
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, createMuonCaloHits(muonCaloHits))

  return pandora::STATUS_CODE_SUCCESS;
}

void DDCaloHitCreatorALLEGRO::getCommonCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                                         PandoraApi::CaloHit::Parameters& caloHitParameters) const {
  const auto position = hit.getPosition();
  const pandora::CartesianVector positionVector(position.x, position.y, position.z);

  // FIXME! AD: for ECAL the cell gemoetry should be pandora::POINTING with cellSize0 = DeltaEta and cellSize1 =
  // DeltaPhi
  caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
  caloHitParameters.m_positionVector = positionVector;
  caloHitParameters.m_expectedDirection = positionVector.GetUnitVector();
  caloHitParameters.m_pParentAddress = static_cast<const void*>(&hit);
  caloHitParameters.m_inputEnergy = hit.getEnergy();
  caloHitParameters.m_time = hit.getTime();
}
