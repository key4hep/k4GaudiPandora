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
 *  @file   k4GaudiPandora/include/GeometryCreatorIdea.h
 *
 *  @brief  Header file for the geometry creator class.
 *
 *  $Log: $
 */

#ifndef GeometryCreatorIdea_h
#define GeometryCreatorIdea_h

#include "Api/PandoraApi.h"

#include "DDGeometryCreator.h"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Geometry creator for the IDEA detector
 */
class GeometryCreatorIdea : public DDGeometryCreator {
public:
  class Settings : public DDGeometryCreator::Settings {
  public:
    Settings() = default;
    ~Settings() = default;

    /// Whether the geometry has a dual-readout HCAL endcap.  The reduced CI geometry is
    /// barrel-only, and registering an endcap subdetector it does not have would fail the
    /// DetType lookup in getExtension.
    bool m_hasHcalEndcap = true;
  };

  /**
   *  @brief  Constructor
   *
   *  @param  settings the creator settings
   *  @param  pPandora address of the relevant pandora instance
   */
  GeometryCreatorIdea(const Settings& settings, pandora::Pandora& pPandora, Gaudi::Algorithm* algorithm);

  /**
   *  @brief  Create geometry
   */
  pandora::StatusCode CreateGeometry() const override;

private:
  /**
   *  @brief  Set mandatory sub detector parameters
   *
   *  @param  subDetectorTypeMap the sub detector type map
   */
  void SetMandatorySubDetectorParameters(SubDetectorTypeMap& subDetectorTypeMap) const override;

  const bool m_hasHcalEndcap; ///< see Settings::m_hasHcalEndcap

  // IDEA ECAL parameters (fiber DRC for o1, crystal DRC for o2)
  void SetEcalParameters(const dd4hep::rec::LayeredCalorimeterData& inputParameters,
                         PandoraApi::Geometry::SubDetector::Parameters& paramBarrel,
                         PandoraApi::Geometry::SubDetector::Parameters& paramEndcap) const;
  void SetHcalBarrelParameters(const dd4hep::rec::LayeredCalorimeterData& inputParameters,
                               PandoraApi::Geometry::SubDetector::Parameters& paramBarrel) const;
  void SetHcalEndcapParameters(const dd4hep::rec::LayeredCalorimeterData& inputParameters,
                               PandoraApi::Geometry::SubDetector::Parameters& paramEndcap) const;
};

#endif
