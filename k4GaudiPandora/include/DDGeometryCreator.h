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

#ifndef DDGEOMETRY_CREATOR_H
#define DDGEOMETRY_CREATOR_H 1

#include <Api/PandoraApi.h>

#include <DDRec/DetectorData.h>

#include <Gaudi/Algorithm.h>

class DDGeometryCreator {
public:
  class Settings {
  public:
    Settings();
    // NN: Material properties obtained from geometry; relevant variables removed
    /// NN: Extra geometry variables removed since all accessible from DDRec
    // NN: Detector names not needed anymore since accessible by type flags
    /// Added by Nikiforos
    bool m_createGaps; ///< Whether should create gaps
  };

  /**
   *  @brief  Constructor
   *
   *  @param  settings the creator settings
   *  @param  pandora reference to the relevant pandora instance
   */
  DDGeometryCreator(const Settings& settings, pandora::Pandora& pandora, Gaudi::Algorithm* algorithm);

  pandora::StatusCode CreateGeometry() const;

protected:
  typedef std::map<pandora::SubDetectorType, PandoraApi::Geometry::SubDetector::Parameters> SubDetectorTypeMap;
  typedef std::map<std::string, PandoraApi::Geometry::SubDetector::Parameters> SubDetectorNameMap;
  /**
   *  @brief  Set mandatory sub detector parameters
   *
   *  @param  subDetectorTypeMap the sub detector type map
   */
  void SetMandatorySubDetectorParameters(SubDetectorTypeMap& subDetectorTypeMap) const;

  /**
   *  @brief  Set additional sub detector parameters
   *
   *  @param  subDetectorNameMap the sub detector name map (for smaller sub detectors, identified uniquely only by name)
   */
  void SetAdditionalSubDetectorParameters(SubDetectorNameMap& subDetectorNameMap) const;

  /**
   *  @brief  Set sub detector parameters to their gear default values
   *
   *  @param  inputParameters input parameters, from gear
   *  @param  subDetectorName the sub detector name
   *  @param  subDetectorType the sub detector type
   *  @param  parameters the sub detector parameters
   */
  void SetDefaultSubDetectorParameters(const dd4hep::rec::LayeredCalorimeterData& inputParameters,
                                       const std::string& subDetectorName,
                                       const pandora::SubDetectorType subDetectorType,
                                       PandoraApi::Geometry::SubDetector::Parameters& parameters) const;

  /**
   *  @brief  Set positions of gaps in ILD detector and add information missing from GEAR parameters file
   *
   *  @param  subDetectorTypeMap the sub detector type map
   *  @param  subDetectorNameMap the sub detector name map (for smaller sub detectors, identified uniquely only by name)
   */
  pandora::StatusCode SetILDSpecificGeometry(SubDetectorTypeMap& subDetectorTypeMap,
                                             SubDetectorNameMap& subDetectorNameMap) const;

  /**
   *  @brief  Add information missing from GEAR parameters file for ILD SDHCAL detector
   *
   *  @param  subDetectorTypeMap the sub detector type map
   */

  /**
   *  @brief  Specify positions of hcal barrel box gaps - ILD specific
   */
  pandora::StatusCode CreateHCalBarrelBoxGaps() const;

  /**
   *  @brief  Specify positions of hcal end cap box gaps - ILD specific
   */
  pandora::StatusCode CreateHCalEndCapBoxGaps() const;

  /**
   *  @brief  Specify positions of hcal barrel concentric polygon gaps - ILD specific
   */
  pandora::StatusCode CreateHCalBarrelConcentricGaps() const;

  /**
   *  @brief  Create box gaps at regular positions on polygonal prism, oriented along main z axis - ILD specific
   *
   *  @param  symmetryOrder the pandora geometry parameters
   *  @param  phi0 the phi coordinate
   *  @param  innerRadius the inner r coordinate
   *  @param  outerRadius the outer r coordinate
   *  @param  minZ the minimum z coordinate
   *  @param  maxZ the maximum z coordinate
   *  @param  gapWidth the gap width
   *  @param  vertexOffset position offset for vertex that doesn't point back to origin of xy plane
   */
  pandora::StatusCode CreateRegularBoxGaps(unsigned int symmetryOrder, float phi0, float innerRadius, float outerRadius,
                                           float minZ, float maxZ, float gapWidth,
                                           pandora::CartesianVector vertexOffset = pandora::CartesianVector(0, 0,
                                                                                                            0)) const;

  const Settings m_settings;    ///< The geometry creator settings
  pandora::Pandora& m_pPandora; ///< Reference to the pandora object to create the geometry
  Gaudi::Algorithm& m_algorithm;
};

#endif // #ifndef GEOMETRY_CREATOR_H
