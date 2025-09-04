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

#ifndef K4GAUDIPANDORA_DDCALO_HIT_CREATOR_H
#define K4GAUDIPANDORA_DDCALO_HIT_CREATOR_H 1

#include <Api/PandoraApi.h>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>

#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/CalorimeterHitCollection.h>

#include <string>
#include <vector>

#include "GaudiKernel/Algorithm.h"

typedef std::vector<edm4hep::CalorimeterHit> CalorimeterHitVector;
typedef std::vector<const edm4hep::CalorimeterHitCollection*> HitCollectionVector;

/**
 *  @brief  DDCaloHitCreator class
 */
class DDCaloHitCreator {
public:
  typedef std::vector<std::string> StringVector;
  typedef std::vector<float> FloatVector;

  /**
   *  @brief  Settings class
   */
  class Settings {
  public:
    /**
     *  @brief  Default constructor
     */
    Settings();

    StringVector m_eCalCaloHitCollections;  ///< The ecal calorimeter hit collections
    StringVector m_hCalCaloHitCollections;  ///< The hcal calorimeter hit collections
    StringVector m_lCalCaloHitCollections;  ///< The lcal calorimeter hit collections
    StringVector m_lHCalCaloHitCollections; ///< The lhcal calorimeter hit collections
    StringVector m_muonCaloHitCollections;  ///< The muon calorimeter hit collections

    float m_eCalToMip;        ///< The calibration from deposited ECal energy to mip
    float m_hCalToMip;        ///< The calibration from deposited HCal energy to mip
    float m_muonToMip;        ///< The calibration from deposited Muon energy to mip
    float m_eCalMipThreshold; ///< Threshold for creating calo hits in the ECal, units mip
    float m_hCalMipThreshold; ///< Threshold for creating calo hits in the HCal, units mip
    float m_muonMipThreshold; ///< Threshold for creating calo hits in the Muon system, units mip

    float m_eCalToEMGeV;        ///< The calibration from deposited ECal energy to EM energy
    float m_eCalToHadGeVBarrel; ///< The calibration from deposited ECal barrel energy to hadronic energy
    float m_eCalToHadGeVEndCap; ///< The calibration from deposited ECal endcap energy to hadronic energy
    float m_hCalToEMGeV;        ///< The calibration from deposited HCal energy to EM energy
    float m_hCalToHadGeV;       ///< The calibration from deposited HCal energy to hadronic energy
    int m_muonDigitalHits;      ///< Muon hits are treated as digital (energy from hit count)
    float m_muonHitEnergy;      ///< The energy for a digital muon calorimeter hit, units GeV

    float m_maxHCalHitHadronicEnergy;      ///< The maximum hadronic energy allowed for a single hcal hit
    int m_nOuterSamplingLayers;            ///< Number of layers from edge for hit to be flagged as an outer layer hit
    float m_layersFromEdgeMaxRearDistance; ///< Maximum number of layers from candidate outer layer hit to rear of
                                           ///< detector

    float m_hCalEndCapInnerSymmetryOrder; ///< HCal end cap inner symmetry order
    float m_hCalEndCapInnerPhiCoordinate; ///< HCal end cap inner phi coordinate

    int m_stripSplittingOn;       ///< Strip splitting activation flag
    int m_useEcalScLayers;        ///< Hybrid ECAL scintillator layers flag
    float m_eCalSiToMip;          ///< Calibration for silicon layers
    float m_eCalScToMip;          ///< Calibration for scintillator layers
    float m_eCalSiMipThreshold;   ///< MIP threshold for silicon layers
    float m_eCalScMipThreshold;   ///< MIP threshold for scintillator layers
    float m_eCalSiToEMGeV;        ///< EM calibration for silicon layers
    float m_eCalScToEMGeV;        ///< EM calibration for scintillator layers
    float m_eCalSiToHadGeVBarrel; ///< Hadronic calibration for silicon layers in the barrel
    float m_eCalScToHadGeVBarrel; ///< Hadronic calibration for scintillator layers in the barrel
    float m_eCalSiToHadGeVEndCap; ///< Hadronic calibration for silicon layers in the endcap
    float m_eCalScToHadGeVEndCap; ///< Hadronic calibration for scintillator layers in the endcap

    float m_eCalBarrelOuterZ; ///< ECal barrel outer Z coordinate
    float m_hCalBarrelOuterZ; ///< HCal barrel outer Z coordinate
    float m_muonBarrelOuterZ; ///< Muon barrel outer Z coordinate
    float m_coilOuterR;       ///< Coil outer radius

    float m_eCalBarrelInnerPhi0;            ///< ECal barrel inner phi0 coordinate
    unsigned int m_eCalBarrelInnerSymmetry; ///< ECal barrel inner symmetry order
    float m_hCalBarrelInnerPhi0;            ///< HCal barrel inner phi0 coordinate
    unsigned int m_hCalBarrelInnerSymmetry; ///< HCal barrel inner symmetry order
    float m_muonBarrelInnerPhi0;            ///< Muon barrel inner phi0 coordinate
    unsigned int m_muonBarrelInnerSymmetry; ///< Muon barrel inner symmetry order

    float m_hCalEndCapOuterR;               ///< HCal endcap outer radius
    float m_hCalEndCapOuterZ;               ///< HCal endcap outer Z coordinate
    float m_hCalBarrelOuterR;               ///< HCal barrel outer radius
    float m_hCalBarrelOuterPhi0;            ///< HCal barrel outer phi0 coordinate
    unsigned int m_hCalBarrelOuterSymmetry; ///< HCal barrel outer symmetry order
    bool m_useSystemId;                     ///< flag whether to use systemId or not to identify origin of the CaloHit
    int m_ecalBarrelSystemId;               ///< systemId of ECal Barrel
    int m_hcalBarrelSystemId;               ///< systemId of HCal Barrel

  public:
    FloatVector m_eCalBarrelNormalVector;
    FloatVector m_hCalBarrelNormalVector;
    FloatVector m_muonBarrelNormalVector;
  };

  /**
   *  @brief  Constructor
   *
   *  @param  pAlgorithm pointer to the Gaudi algorithm instance
   *  @param  settings the creator settings
   *  @param  pandora reference to the relevant pandora instance
   */
  DDCaloHitCreator(const Settings& settings, pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm);
  virtual ~DDCaloHitCreator();

  pandora::StatusCode
  createECalCaloHits(const std::map<std::string, std::vector<edm4hep::CalorimeterHit>>& inputECalCaloHits) const;
  pandora::StatusCode createHCalCaloHits(const std::vector<edm4hep::CalorimeterHit>& hCalCaloHits) const;
  pandora::StatusCode createMuonCaloHits(const std::vector<edm4hep::CalorimeterHit>& muonCaloHits) const;
  pandora::StatusCode createLCalCaloHits(const std::vector<edm4hep::CalorimeterHit>& lCalCaloHits) const;
  pandora::StatusCode createLHCalCaloHits(const std::vector<edm4hep::CalorimeterHit>& LHCalCaloHits) const;

  /**
   *  @brief  Create calo hits
   *
   *  @param  inputHits
   *  @param  outputHits
   */
  virtual pandora::StatusCode
  createCaloHits(const std::map<std::string, std::vector<edm4hep::CalorimeterHit>>& eCaloHitsMap,
                 const std::vector<edm4hep::CalorimeterHit>& hCalCaloHits,
                 const std::vector<edm4hep::CalorimeterHit>& muonCaloHits,
                 const std::vector<edm4hep::CalorimeterHit>& lCalCaloHits,
                 const std::vector<edm4hep::CalorimeterHit>& lhCalCaloHits) const;

  /**
   *  @brief  Get the calorimeter hit vector
   *
   *  @return The calorimeter hit vector
   */
  const CalorimeterHitVector& GetCalorimeterHitVector() const;

  /**
   *  @brief  Reset the calo hit creator
   */
  void Reset();

protected:
  /**
   *  @brief  Get common calo hit properties: position, parent address, input energy and time
   *
   *  @param  pCaloHit the lcio calorimeter hit
   *  @param  caloHitParameters the calo hit parameters to populate
   */
  void getCommonCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                  PandoraApi::CaloHit::Parameters& caloHitParameters) const;

  /**
   *  @brief  Get end cap specific calo hit properties: cell size, absorber radiation and interaction lengths, normal
   * vector
   *
   *  @param  pCaloHit the lcio calorimeter hit
   *  @param  layers the vector of layers from DDRec extensions
   *  @param  caloHitParameters the calo hit parameters to populate
   *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
   */
  void getEndCapCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
                                  PandoraApi::CaloHit::Parameters& caloHitParameters, float& absorberCorrection) const;

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
  void getBarrelCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
                                  unsigned int barrelSymmetryOrder, PandoraApi::CaloHit::Parameters& caloHitParameters,
                                  const std::vector<float>& normalVector, float& absorberCorrection) const;

  /**
   *  @brief  Get number of active layers from position of a calo hit to the edge of the detector
   *
   *  @param  pCaloHit the lcio calorimeter hit
   */
  int getNLayersFromEdge(const edm4hep::CalorimeterHit& hit) const;

  /**
   *  @brief  Get the maximum radius of a calo hit in a polygonal detector structure
   *
   *  @param  pCaloHit the lcio calorimeter hit
   *  @param  symmetryOrder the symmetry order
   *  @param  phi0 the angular orientation
   *
   *  @return the maximum radius
   */
  float getMaximumRadius(const edm4hep::CalorimeterHit& pCaloHit, const unsigned int symmetryOrder,
                         const float phi0) const;

  const Settings m_settings; ///< The calo hit creator settings

  pandora::Pandora& m_pandora; ///< Reference to the pandora object to create calo hits

  float m_hCalBarrelLayerThickness; ///< HCal barrel layer thickness
  float m_hCalEndCapLayerThickness; ///< HCal endcap layer thickness

  CalorimeterHitVector m_calorimeterHitVector; ///< The calorimeter hit vector

  dd4hep::VolumeManager m_volumeManager; ///< DD4hep volume manager

  const Gaudi::Algorithm& m_algorithm; ///< Pointer to the Gaudi algorithm for logging
};

#endif // K4GAUDIPANDORA_DDCALO_HIT_CREATOR_H
