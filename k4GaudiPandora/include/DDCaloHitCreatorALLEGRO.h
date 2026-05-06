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
 *  @file   DDMarlinPandora/include/DDCaloHitCreatorALLEGRO.h
 *
 *  @brief  Header file for the calo hit creator class.
 *
 *  $Log: $
 */

#ifndef DDCALO_HIT_CREATORALLEGRO_H
#define DDCALO_HIT_CREATORALLEGRO_H 1

#include "Api/PandoraApi.h"

#include "DDCaloHitCreator.h"

/**
 *  @brief  DDCaloHitCreator class
 */
class DDCaloHitCreatorALLEGRO : public DDCaloHitCreator {
public:
  /**
   *  @brief  Constructor
   *
   *  @param  settings the creator settings
   *  @param  pPandora address of the relevant pandora instance
   */
  DDCaloHitCreatorALLEGRO(const Settings& settings, pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm);
  ~DDCaloHitCreatorALLEGRO() override = default;

  /**
   *  @brief  Create calo hits
   *
   *  @param  pLCEvent the lcio event
   */
  virtual pandora::StatusCode createCaloHits(const std::vector<CollectionDescriptor>& eCalCaloHits,
                                             const std::vector<CollectionDescriptor>& hCalCaloHits,
                                             const std::vector<CollectionDescriptor>& muonCaloHits,
                                             const std::vector<edm4hep::CalorimeterHit>& lCalCaloHits,
                                             const std::vector<edm4hep::CalorimeterHit>& lhCalCaloHits) const override;

private:
  /**
   *  @brief  Get common calo hit properties: position, parent address, input energy and time
   *
   *  @param  pCaloHit the lcio calorimeter hit
   *  @param  caloHitParameters the calo hit parameters to populate
   */
  void getCommonCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                  PandoraApi::CaloHit::Parameters& caloHitParameters) const override;

  /**
   *  @brief  Get end cap specific calo hit properties: cell size, absorber radiation and interaction lengths, normal
   * vector
   *
   *  @param  pCaloHit the lcio calorimeter hit
   *  @param  layers the vector of layers from DDRec extensions
   *  @param  caloHitParameters the calo hit parameters to populate
   *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
   */
  virtual void getEndCapCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                          const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
                                          PandoraApi::CaloHit::Parameters& caloHitParameters,
                                          float& absorberCorrection) const override;

  /**
   *  @brief  Get barrel specific calo hit properties: cell size, absorber radiation and interaction lengths, normal
   * vector
   *
   *  @param  pCaloHit the lcio calorimeter hit
   *  @param  layers the vector of layers from DDRec extensions
   *  @param  barrelSymmetryOrder the barrel order of symmetry
   *  @param  caloHitParameters the calo hit parameters to populate
   *  @param  normalVector is the normalVector to the sensitive layers in local coordinates
   *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
   */
  virtual void getBarrelCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                          const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
                                          unsigned int barrelSymmetryOrder,
                                          PandoraApi::CaloHit::Parameters& caloHitParameters,
                                          const std::vector<float>& normalVector,
                                          float& absorberCorrection) const override;

  std::shared_ptr<dd4hep::DDSegmentation::Segmentation> m_hcalBarrelSegmentation = {};
  std::shared_ptr<dd4hep::DDSegmentation::Segmentation> m_hcalEndcapSegmentation = {};
};
#endif
