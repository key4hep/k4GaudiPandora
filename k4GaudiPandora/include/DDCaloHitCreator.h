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
 *  @file   DDMarlinPandora/include/DDCaloHitCreator.h
 *
 *  @brief  Header file for the calo hit creator class.
 *
 *  $Log: $
 */

#ifndef DDCALO_HIT_CREATOR_H
#define DDCALO_HIT_CREATOR_H 1

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHit.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "Api/PandoraApi.h"
#include <vector>
#include <string>

/**
 *  @brief  DDCaloHitCreator class
 */
class DDCaloHitCreator {
public:
  typedef std::vector<std::string> StringVector;
  typedef std::vector<float>       FloatVector;

  /**
   *  @brief  Settings class
   */
  class Settings {
  public:
    /**
     *  @brief  Default constructor
     */
    Settings();

    StringVector m_eCalCaloHitCollections;   ///< The ecal calorimeter hit collections
    StringVector m_hCalCaloHitCollections;   ///< The hcal calorimeter hit collections
    StringVector m_lCalCaloHitCollections;   ///< The lcal calorimeter hit collections
    StringVector m_lHCalCaloHitCollections;  ///< The lhcal calorimeter hit collections
    StringVector m_muonCaloHitCollections;   ///< The muon calorimeter hit collections

    float m_eCalToMip;         ///< The calibration from deposited ECal energy to mip
    float m_hCalToMip;         ///< The calibration from deposited HCal energy to mip
    float m_muonToMip;         ///< The calibration from deposited Muon energy to mip
    float m_eCalMipThreshold;  ///< Threshold for creating calo hits in the ECal, units mip
    float m_hCalMipThreshold;  ///< Threshold for creating calo hits in the HCal, units mip
    float m_muonMipThreshold;  ///< Threshold for creating calo hits in the Muon system, units mip

    float m_eCalToEMGeV;         ///< The calibration from deposited ECal energy to EM energy
    float m_eCalToHadGeVBarrel;  ///< The calibration from deposited ECal barrel energy to hadronic energy
    float m_eCalToHadGeVEndCap;  ///< The calibration from deposited ECal endcap energy to hadronic energy
    float m_hCalToEMGeV;         ///< The calibration from deposited HCal energy to EM energy
    float m_hCalToHadGeV;        ///< The calibration from deposited HCal energy to hadronic energy
    int   m_muonDigitalHits;     ///< Muon hits are treated as digital (energy from hit count)
    float m_muonHitEnergy;       ///< The energy for a digital muon calorimeter hit, units GeV

    float m_maxHCalHitHadronicEnergy;  ///< The maximum hadronic energy allowed for a single hcal hit
    int   m_nOuterSamplingLayers;      ///< Number of layers from edge for hit to be flagged as an outer layer hit
    float m_layersFromEdgeMaxRearDistance;  ///< Maximum number of layers from candidate outer layer hit to rear of detector

    float m_hCalEndCapInnerSymmetryOrder;  ///< HCal end cap inner symmetry order
    float m_hCalEndCapInnerPhiCoordinate;  ///< HCal end cap inner phi coordinate

    int   m_stripSplittingOn;      ///< Strip splitting activation flag
    int   m_useEcalScLayers;       ///< Hybrid ECAL scintillator layers flag
    float m_eCalSiToMip;           ///< Calibration for silicon layers
    float m_eCalScToMip;           ///< Calibration for scintillator layers
    float m_eCalSiMipThreshold;    ///< MIP threshold for silicon layers
    float m_eCalScMipThreshold;    ///< MIP threshold for scintillator layers
    float m_eCalSiToEMGeV;         ///< EM calibration for silicon layers
    float m_eCalScToEMGeV;         ///< EM calibration for scintillator layers
    float m_eCalSiToHadGeVBarrel;  ///< Hadronic calibration for silicon layers in the barrel
    float m_eCalScToHadGeVBarrel;  ///< Hadronic calibration for scintillator layers in the barrel
    float m_eCalSiToHadGeVEndCap;  ///< Hadronic calibration for silicon layers in the endcap
    float m_eCalScToHadGeVEndCap;  ///< Hadronic calibration for scintillator layers in the endcap

    float m_eCalBarrelOuterZ;  ///< ECal barrel outer Z coordinate
    float m_hCalBarrelOuterZ;  ///< HCal barrel outer Z coordinate
    float m_muonBarrelOuterZ;  ///< Muon barrel outer Z coordinate
    float m_coilOuterR;        ///< Coil outer radius

    float        m_eCalBarrelInnerPhi0;      ///< ECal barrel inner phi0 coordinate
    unsigned int m_eCalBarrelInnerSymmetry;  ///< ECal barrel inner symmetry order
    float        m_hCalBarrelInnerPhi0;      ///< HCal barrel inner phi0 coordinate
    unsigned int m_hCalBarrelInnerSymmetry;  ///< HCal barrel inner symmetry order
    float        m_muonBarrelInnerPhi0;      ///< Muon barrel inner phi0 coordinate
    unsigned int m_muonBarrelInnerSymmetry;  ///< Muon barrel inner symmetry order

    float        m_hCalEndCapOuterR;         ///< HCal endcap outer radius
    float        m_hCalEndCapOuterZ;         ///< HCal endcap outer Z coordinate
    float        m_hCalBarrelOuterR;         ///< HCal barrel outer radius
    float        m_hCalBarrelOuterPhi0;      ///< HCal barrel outer phi0 coordinate
    unsigned int m_hCalBarrelOuterSymmetry;  ///< HCal barrel outer symmetry order
  };

  DDCaloHitCreator(const Settings& settings, const pandora::Pandora* const pPandora);
  ~DDCaloHitCreator();

  pandora::StatusCode CreateECalCaloHits(const edm4hep::CalorimeterHitCollection& inputECalCaloHits);
  pandora::StatusCode CreateHCalCaloHits(const edm4hep::CalorimeterHitCollection& inputHCalCaloHits);
  pandora::StatusCode CreateMuonCaloHits(const edm4hep::CalorimeterHitCollection& inputMuonCaloHits);
  pandora::StatusCode CreateLCalCaloHits(const edm4hep::CalorimeterHitCollection& inputLCalCaloHits);
  pandora::StatusCode CreateLHCalCaloHits(const edm4hep::CalorimeterHitCollection& inputLHCalCaloHits);

  const pandora::Pandora& m_pandora;
};

#endif  // DDCALO_HIT_CREATOR_H

