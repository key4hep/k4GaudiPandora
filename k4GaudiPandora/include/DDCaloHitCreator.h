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

#include <edm4hep/CalorimeterHit.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/MutableCalorimeterHit.h>

#include <k4FWCore/Transformer.h>

#include "DDSegmentation/BitFieldCoder.h"

#include "Api/PandoraApi.h"

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>

typedef std::vector<std::shared_ptr<edm4hep::CalorimeterHit>> CalorimeterHitVector;

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
    //NN: Material properties variables removed; obtained from geometry

    float m_eCalToMip;         ///< The calibration from deposited ECal energy to mip
    float m_hCalToMip;         ///< The calibration from deposited HCal energy to mip
    float m_muonToMip;         ///< The calibration from deposited Muon energy to mip
    float m_eCalMipThreshold;  ///< Threshold for creating calo hits in the ECal, units mip
    float m_hCalMipThreshold;  ///< Threshold for creating calo hits in the HCal, units mip
    float m_muonMipThreshold;  ///< Threshold for creating calo hits in the HCal, units mip

    float m_eCalToEMGeV;         ///< The calibration from deposited ECal energy to EM energy
    float m_eCalToHadGeVBarrel;  ///< The calibration from deposited ECal barrel energy to hadronic energy
    float m_eCalToHadGeVEndCap;  ///< The calibration from deposited ECal endcap energy to hadronic energy
    float m_hCalToEMGeV;         ///< The calibration from deposited HCal energy to EM energy
    float m_hCalToHadGeV;        ///< The calibration from deposited HCal energy to hadronic energy
    int   m_muonDigitalHits;     ///< Muon hits are treated as digital (energy from hit count)
    float m_muonHitEnergy;       ///< The energy for a digital muon calorimeter hit, units GeV

    float m_maxHCalHitHadronicEnergy;  ///< The maximum hadronic energy allowed for a single hcal hit
    int   m_nOuterSamplingLayers;      ///< Number of layers from edge for hit to be flagged as an outer layer hit
    float
        m_layersFromEdgeMaxRearDistance;  ///< Maximum number of layers from candidate outer layer hit to rear of detector

    int   m_hCalEndCapInnerSymmetryOrder;  ///< HCal end cap inner symmetry order (missing from ILD00 gear file)
    float m_hCalEndCapInnerPhiCoordinate;  ///< HCal end cap inner phi coordinate (missing from ILD00 gear file)

    // For Strip Splitting method and hybrid ECAL.
    int   m_stripSplittingOn;      ///< To use SSA, this should be true (default is false)
    std::vector<int> m_useEcalScLayers;       ///< To use scintillator layers ~ hybrid ECAL, this should be true (default is false)
    std::vector<int> m_useEcalSiLayers;       ///< To use scintillator layers ~ hybrid ECAL, this should be true (default is false)
    float m_eCalSiToMip;           ///< The calibration from deposited Si-layer energy to mip
    float m_eCalScToMip;           ///< The calibration from deposited Sc-layer energy to mip
    float m_eCalSiMipThreshold;    ///< Threshold for creating calo hits in the Si-layers of ECAL, units mip
    float m_eCalScMipThreshold;    ///< Threshold for creating calo hits in the Sc-layers of ECAL, units mip
    float m_eCalSiToEMGeV;         ///< The calibration from deposited Si-layer energy to EM energy
    float m_eCalScToEMGeV;         ///< The calibration from deposited Sc-layer energy to EM energy
    float m_eCalSiToHadGeVBarrel;  ///< The calibration from deposited Si-layer energy on the enecaps to hadronic energy
    float m_eCalScToHadGeVBarrel;  ///< The calibration from deposited Sc-layer energy on the endcaps to hadronic energy
    float m_eCalSiToHadGeVEndCap;  ///< The calibration from deposited Si-layer energy on the enecaps to hadronic energy
    float m_eCalScToHadGeVEndCap;  ///< The calibration from deposited Sc-layer energy on the endcaps to hadronic energy

    ///ADDED BY NIKIFOROS
    //Note that names are not needed anymore since the detector elements can be accessed by type flags
    //Nikiforos: Moved from main class

    float m_eCalBarrelOuterZ;  ///< ECal barrel outer z coordinate
    float m_hCalBarrelOuterZ;  ///< HCal barrel outer z coordinate
    float m_muonBarrelOuterZ;  ///< Muon barrel outer z coordinate
    float m_coilOuterR;        ///< Coil outer r coordinate

    float        m_eCalBarrelInnerPhi0;      ///< ECal barrel inner phi0 coordinate
    unsigned int m_eCalBarrelInnerSymmetry;  ///< ECal barrel inner symmetry order
    float        m_hCalBarrelInnerPhi0;      ///< HCal barrel inner phi0 coordinate
    unsigned int m_hCalBarrelInnerSymmetry;  ///< HCal barrel inner symmetry order
    float        m_muonBarrelInnerPhi0;      ///< Muon barrel inner phi0 coordinate
    unsigned int m_muonBarrelInnerSymmetry;  ///< Muon barrel inner symmetry order

    float        m_hCalEndCapOuterR;         ///< HCal endcap outer r coordinate
    float        m_hCalEndCapOuterZ;         ///< HCal endcap outer z coordinate
    float        m_hCalBarrelOuterR;         ///< HCal barrel outer r coordinate
    float        m_hCalBarrelOuterPhi0;      ///< HCal barrel outer phi0 coordinate
    unsigned int m_hCalBarrelOuterSymmetry;  ///< HCal barrel outer symmetry order

    std::string m_caloEncodingString; ///< Encoding String for calo

  public:
    FloatVector m_eCalBarrelNormalVector;
    FloatVector m_hCalBarrelNormalVector;
    FloatVector m_muonBarrelNormalVector;
  };

  /**
     *  @brief  Constructor
     *
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
  DDCaloHitCreator(const Settings& settings, const pandora::Pandora* const pPandora, IMessageSvc* msgSvc);

  /**
     *  @brief  Destructor
     */
  ~DDCaloHitCreator();

  /**
     *  @brief  Create calo hits
     *
     *  @param  eCalCollections  CalorimeterHit Collection for the ECal
     *  @param  hCalCollections  CalorimeterHit Collection for the HCal
     *  @param  mCalCollections  CalorimeterHit Collection for the Muon Calo
     *  @param  lCalCollections  CalorimeterHit Collection for the LCal
     *  @param  lhCalCollections CalorimeterHit Collection for the HLCal
     *  @return The calorimeter hit vector
     */
  pandora::StatusCode CreateCaloHits(
   const std::vector<const edm4hep::CalorimeterHitCollection*>& eCalCollections,
   const std::vector<const edm4hep::CalorimeterHitCollection*>& hCalCollections,
   const std::vector<const edm4hep::CalorimeterHitCollection*>& mCalCollections,
   const std::vector<const edm4hep::CalorimeterHitCollection*>& lCalCollections,
   const std::vector<const edm4hep::CalorimeterHitCollection*>& lhCalCollections
 );

  /**
     *  @brief  Get the calorimeter hit vector
     */
  const CalorimeterHitVector& GetCalorimeterHitVector() const;

  /**
     *  @brief  Reset the calo hit creator
     */
  void Reset();

private:
  /**
     *  @brief  Create ecal calo hits
     *
     *  @param  eCalCollections the eCal collections
     */
  pandora::StatusCode CreateECalCaloHits(const std::vector<const edm4hep::CalorimeterHitCollection*>& eCalCollections);

  /**
     *  @brief  Create hcal calo hits
     *
     *  @param  hCalCollections the hCal collections
     */
  pandora::StatusCode CreateHCalCaloHits(const std::vector<const edm4hep::CalorimeterHitCollection*>& hCalCollections);

  /**
     *  @brief  Create muon calo hits
     *
     *  @param  mCalCollections the mCal collections
     */
  pandora::StatusCode CreateMuonCaloHits(const std::vector<const edm4hep::CalorimeterHitCollection*>& mCalCollections);

  /**
     *  @brief  Create lcal calo hits
     *
     *  @param  lCalCollections the lCal collections
     */
  pandora::StatusCode CreateLCalCaloHits(const std::vector<const edm4hep::CalorimeterHitCollection*>& mCalCollections);

  /**
     *  @brief  Create lhcal calo hits
     *
     *  @param  lhCalCollections the lhCal collections
     */
  pandora::StatusCode CreateLHCalCaloHits(const std::vector<const edm4hep::CalorimeterHitCollection*>& lhCalCollections);

  /**
     *  @brief  Get common calo hit properties: position, parent address, input energy and time
     *
     *  @param  pCaloHit the lcio calorimeter hit
     *  @param  caloHitParameters the calo hit parameters to populate
     */
  void GetCommonCaloHitProperties(std::shared_ptr<edm4hep::CalorimeterHit>,
                                  PandoraApi::CaloHit::Parameters&   caloHitParameters) const;

  /**
     *  @brief  Get end cap specific calo hit properties: cell size, absorber radiation and interaction lengths, normal vector
     *
     *  @param  pCaloHit the edm4hep calorimeter hit
     *  @param  layers the vector of layers from DDRec extensions
     *  @param  caloHitParameters the calo hit parameters to populate
     *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
     */
  void GetEndCapCaloHitProperties(std::shared_ptr<edm4hep::CalorimeterHit> pCaloHit,
                                  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
                                  PandoraApi::CaloHit::Parameters& caloHitParameters, float& absorberCorrection) const;

  /**
     *  @brief  Get barrel specific calo hit properties: cell size, absorber radiation and interaction lengths, normal vector
     *
     *  @param  pCaloHit the edm4hep calorimeter hit
     *  @param  layers the vector of layers from DDRec extensions
     *  @param  barrelSymmetryOrder the barrel order of symmetry
     *  @param  caloHitParameters the calo hit parameters to populate
     *  @param  normalVector is the normalVector to the sensitive layers in local coordinates
     *  @param  absorberCorrection to receive the absorber thickness correction for the mip equivalent energy
     */
  void GetBarrelCaloHitProperties(std::shared_ptr<edm4hep::CalorimeterHit> pCaloHit,
                                  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& layers,
                                  unsigned int barrelSymmetryOrder, PandoraApi::CaloHit::Parameters& caloHitParameters,
                                  FloatVector const& normalVector, float& absorberCorrection) const;

  /**
     *  @brief  Get number of active layers from position of a calo hit to the edge of the detector
     *
     *  @param  pCaloHit the edm4hep calorimeter hit
     */
  int GetNLayersFromEdge(std::shared_ptr<edm4hep::CalorimeterHit> pCaloHit) const;

  /**
     *  @brief  Get the maximum radius of a calo hit in a polygonal detector structure
     *
     *  @param  pCaloHit the edm4hep calorimeter hit
     *  @param  symmetryOrder the symmetry order
     *  @param  phi0 the angular orientation
     *
     *  @return the maximum radius
     */
  float GetMaximumRadius(std::shared_ptr<edm4hep::CalorimeterHit> pCaloHit, const unsigned int symmetryOrder,
                         const float phi0) const;

  const Settings m_settings;  ///< The calo hit creator settings

  const pandora::Pandora& m_pandora;  ///< Reference to the pandora object to create calo hits

  float m_hCalBarrelLayerThickness;  ///< HCal barrel layer thickness
  float m_hCalEndCapLayerThickness;  ///< HCal endcap layer thickness

  CalorimeterHitVector m_calorimeterHitVector;  ///< The calorimeter hit vector

  dd4hep::VolumeManager m_volumeManager;  ///< DD4hep volume manager

  dd4hep::DDSegmentation::BitFieldCoder m_cell_encoder;
  IMessageSvc* m_msgSvc;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const CalorimeterHitVector& DDCaloHitCreator::GetCalorimeterHitVector() const { return m_calorimeterHitVector; }

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DDCaloHitCreator::Reset() { m_calorimeterHitVector.clear(); }

#endif  // #ifndef CALO_HIT_CREATOR_H
