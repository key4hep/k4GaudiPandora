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
#ifndef DDCALODIGI_H
#define DDCALODIGI_H 1

#include "DDRec/DetectorData.h"
#include "Gaudi/Accumulators/RootHistogram.h"
#include "Gaudi/Property.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "CalorimeterHitType.h"
#include "DDScintillatorPpdDigi.h"

#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RandGauss.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <random>
#include <string>
#include <vector>

/** === DDCaloDigi Processor === <br>
 *  Simple calorimeter digitizer Processor. <br>
 *  Ported from ILDCaloDigi to use DD4hep
 *  Takes SimCalorimeterHit Collections and <br>
 *  produces CalorimeterHit Collections. <br>
 *  Simulated energy depositions in active <br>
 *  layers of calorimeters are <br>
 *  converted into physical energy. This is done <br>
 *  taking into account sampling fractions of <br>
 *  ECAL and HCAL. <br>
 *  User has to specify ECAL and HCAL SimCalorimeterHit <br>
 *  collections with processor parameters <br>
 *  HCALCollections and ECALCollections. <br>
 *  The names of the output CalorimeterHit Collections <br>
 *  are specified with processor parameters <br>
 *  ECALOutputCollection and HCALOutputCollection. <br>
 *  Conversion factors for ECAL and HCAL <br>
 *  are specified via processor parameters  <br>
 *  CalibrECAL and CalibrHCAL. <br>
 *  It should be noted that ECAL and HCAL may consist <br>
 *  of several sections with different sampling fractions. <br>
 *  To handle this situation, calibration coefficients for <br>
 *  ECAL and HCAL are passed as arrays of floats with each element <br>
 *  in this array corresponding to certain section with <br>
 *  a given sampling fraction. <br>
 *  List of layer numbers terminating each section are given through <br>
 *  processor parameters ECALLayers and HCALLayers <br>
 *  There is an option to perform digitization of <br>
 *  both ECAL and HCAL in a digital mode. <br>
 *  Digital mode is activated by  <br>
 *  setting processor parameters <br>
 *  IfDigitalEcal / IfDigitalHcal to 1. <br>
 *  In this case CalibrECAL / CalibrHCAL will  <br>
 *  convert the number of hits into physical energy. <br>
 *  Thresholds on hit energies in ECAL and HCAL <br>
 *  are set with processor parameters <br>
 *  ECALThreshold and HCALThreshold.  <br>
 *  Relations between CalorimeterHits and SimCalorimeterHits <br>
 *  are held in the corresponding relation collection. <br>
 *  The name of this relation collection is specified <br>
 *  via processor parameter RelationOutputCollection. <br>
 *  <h4>Input collections and prerequisites</h4>
 *  SimCalorimeterHit collections <br>
 *  <h4>Output</h4>
 *  CalorimeterHit collections for ECal and HCal. <br>
 *  Collection of relations <br>
 *  between CalorimeterHits and SimCalorimeterHits. <br>
 *  For ECal Calorimeter hits the variable type is set to 0, <br>
 *  whereas for HCal Calorimeter hits the type is set to 1 <br>
 *  @author A. Raspereza (DESY) <br>
 *  @author M. Thomson (DESY) <br>
 *  @version $Id$ <br>
 */

const int MAX_LAYERS = 200;
const int MAX_STAVES = 16;

using retType = std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>;

using DDCaloDigi_t = k4FWCore::MultiTransformer<retType(const edm4hep::SimCalorimeterHitCollection&,
                                                        const edm4hep::EventHeaderCollection&)>;

struct DDCaloDigi final : DDCaloDigi_t {
  DDCaloDigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  retType operator()(const edm4hep::SimCalorimeterHitCollection& simCaloHits,
                     const edm4hep::EventHeaderCollection&       headers) const;

  //APS: AFAIU we need a registerCallBack function to make the creation of histograms possible
  void registerCallBack(Gaudi::StateMachine::Transition, std::function<void()>) {}

private:
  // check if input collection is ECAL or HCAL
  Gaudi::Property<bool> m_inputColIsECAL{
      this, "InputColIsECAL", {true}, "Input SimCalorimeterHit collection is ECAL? If false, then is HCAL."};

  // digitazing parameters for ECAL and HCAL
  Gaudi::Property<float>       m_thresholdEcal{this, "ECALThreshold", {5.0e-5}, "Threshold for ECAL Hits in GeV"};
  Gaudi::Property<std::string> m_unitThresholdEcal{
      this,
      "ECALThresholdUnit",
      {"GeV"},
      "Unit for ECAL Threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  Gaudi::Property<std::vector<float>> m_thresholdHcal{
      this, "HCALThreshold", {0.00004}, "Threshold for ECAL Hits in GeV"};
  Gaudi::Property<std::string> m_unitThresholdHcal{
      this,
      "HCALThresholdUnit",
      {"GeV"},
      "Unit for HCAL Threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  Gaudi::Property<std::vector<int>>   m_ecalLayers{this, "ECALLayers", {20, 100}, "Index of ECAL Layers"};
  Gaudi::Property<std::vector<int>>   m_hcalLayers{this, "HCALLayers", {100}, "Index of HCAL Layers"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffEcal{
      this, "CalibrECAL", {40.91, 81.81}, "Calibration coefficients for ECAL"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffHcalBarrel{
      this, "CalibrHCALBarrel", {0.0}, "Calibration coefficients for Barrel HCAL"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffHcalEndcap{
      this, "CalibrHCALEndcap", {0.0}, "Calibration coefficients for Endcap HCAL"};
  Gaudi::Property<std::vector<float>> m_calibrCoeffHcalOther{
      this, "CalibrHCALOther", {0.0}, "Calibration coefficients for Other HCAL"};
  Gaudi::Property<int> m_digitalEcal{this, "IfDigitalEcal", {0}, "Digital ECAL"};
  Gaudi::Property<int> m_mapsEcalCorrection{
      this, "MapsEcalCorrection", {0}, "Ecal correction for theta dependency of calibration for MAPS"};
  Gaudi::Property<int>   m_digitalHcal{this, "IfDigitalHcal", {0}, "Digital Hcal"};
  Gaudi::Property<int>   m_ecalGapCorrection{this, "ECALGapCorrection", {1}, "Correct for ECAL gaps"};
  Gaudi::Property<int>   m_hcalGapCorrection{this, "HCALGapCorrection", {1}, "Correct for HCAL gaps"};
  Gaudi::Property<float> m_ecalEndcapCorrectionFactor{
      this, "ECALEndcapCorrectionFactor", {1.025}, "Energy correction for ECAL endcap"};
  Gaudi::Property<float> m_hcalEndcapCorrectionFactor{
      this, "HCALEndcapCorrectionFactor", {1.025}, "Energy correction for HCAL endcap"};
  Gaudi::Property<float> m_ecalGapCorrectionFactor{
      this, "ECALGapCorrectionFactor", {1.0}, "Factor applied to gap correction"};
  Gaudi::Property<float> m_ecalModuleGapCorrectionFactor{
      this, "ECALModuleGapCorrectionFactor", {0.5}, "Factor applied to module gap correction ECAL"};
  Gaudi::Property<float> m_hcalModuleGapCorrectionFactor{
      this, "HCALModuleGapCorrectionFactor", {0.5}, "Factor applied to module gap correction HCAL"};

  // histograms that will be filled --- <1> or <2> is the number of dimensions of the histogram (1D or 2D)
  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fEcal{this, "fEcal", "Ecal time profile", {1000, 0., 1000.0}};
  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fHcal{this, "fHcal", "Hcal time profile", {1000, 0., 1000.0}};

  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fEcalC{
      this, "fEcalC", "Ecal time profile cor", {1000, 0., 1000.0}};
  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fHcalC{
      this, "fHcalC", "Hcal time profile cor", {1000, 0., 1000.0}};

  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fEcalC1{
      this, "fEcalC1", "Ecal time profile cor", {100, 0., 1000.0}};
  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fHcalC1{
      this, "fHcalC1", "Hcal time profile cor", {100, 0., 1000.0}};

  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fEcalC2{
      this, "fEcalC2", "Ecal time profile cor", {10, 0., 1000.0}};
  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fHcalC2{
      this, "fHcalC2", "Hcal time profile cor", {10, 0., 1000.0}};

  mutable Gaudi::Accumulators::StaticHistogram<2> fHcalCvsE{
      this, "fHcalCvsE", "Hcal time profile cor", {{100, 0., 500.0}, {100, 0., 10.}}};

  mutable Gaudi::Accumulators::StaticHistogram<2> fHcalLayer1{
      this, "fHcalLayer1", "Hcal layer 1 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::StaticHistogram<2> fHcalLayer11{
      this, "fHcalLayer11", "Hcal layer 11 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::StaticHistogram<2> fHcalLayer21{
      this, "fHcalLayer21", "Hcal layer 21 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::StaticHistogram<2> fHcalLayer31{
      this, "fHcalLayer31", "Hcal layer 31 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::StaticHistogram<2> fHcalLayer41{
      this, "fHcalLayer41", "Hcal layer 41 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::StaticHistogram<2> fHcalLayer51{
      this, "fHcalLayer51", "Hcal layer 51 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::StaticHistogram<2> fHcalLayer61{
      this, "fHcalLayer61", "Hcal layer 61 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};
  mutable Gaudi::Accumulators::StaticHistogram<2> fHcalLayer71{
      this, "fHcalLayer71", "Hcal layer 71 map", {{300, -4500., 4500.0}, {300, -4500, 4500.}}};

  mutable Gaudi::Accumulators::StaticHistogram<1> fHcalRLayer1{
      this, "fHcalRLayer1", "Hcal R layer 1", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticHistogram<1> fHcalRLayer11{
      this, "fHcalRLayer11", "Hcal R layer 11", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticHistogram<1> fHcalRLayer21{
      this, "fHcalRLayer21", "Hcal R layer 21", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticHistogram<1> fHcalRLayer31{
      this, "fHcalRLayer31", "Hcal R layer 31", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticHistogram<1> fHcalRLayer41{
      this, "fHcalRLayer41", "Hcal R layer 41", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticHistogram<1> fHcalRLayer51{
      this, "fHcalRLayer51", "Hcal R layer 51", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticHistogram<1> fHcalRLayer61{
      this, "fHcalRLayer61", "Hcal R layer 61", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticHistogram<1> fHcalRLayer71{
      this, "fHcalRLayer71", "Hcal R layer 71", {50, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fHcalRLayerNorm{
      this, "fHcalRLayerNorm", "Hcal R layer Norm", {50, 0., 5000.0}};

  mutable Gaudi::Accumulators::StaticHistogram<2> fEcalLayer1{
      this, "fEcalLayer1", "Ecal layer 1 map", {{1800, -4500., 4500.0}, {1800, -4500, 4500.}}};
  mutable Gaudi::Accumulators::StaticHistogram<2> fEcalLayer11{
      this, "fEcalLayer11", "Ecal layer 11 map", {{1800, -4500., 4500.0}, {1800, -4500, 4500.}}};
  mutable Gaudi::Accumulators::StaticHistogram<2> fEcalLayer21{
      this, "fEcalLayer21", "Ecal layer 21 map", {{1800, -4500., 4500.0}, {1800, -4500, 4500.}}};

  mutable Gaudi::Accumulators::StaticHistogram<1> fEcalRLayer1{
      this, "fEcalRLayer1", "Ecal R layer 1", {100, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticHistogram<1> fEcalRLayer11{
      this, "fEcalRLayer11", "Ecal R layer 11", {100, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticHistogram<1> fEcalRLayer21{
      this, "fEcalRLayer21", "Ecal R layer 21", {100, 0., 5000.0}};
  mutable Gaudi::Accumulators::StaticWeightedHistogram<1> fEcalRLayerNorm{
      this, "fEcalRLayerNorm", "Ecal R layer Norm", {100, 0., 5000.0}};

  // timing parameters for ECAL
  Gaudi::Property<int> m_useEcalTiming{this, "UseEcalTiming", {0}, "Use ECAL hit times"};
  Gaudi::Property<int> m_ecalCorrectTimesForPropagation{
      this, "ECALCorrectTimesForPropagation", {0}, "Correct ECAL hit times for propagation: radial distance/c"};
  Gaudi::Property<float> m_ecalTimeWindowMin{this, "ECALTimeWindowMin", {-10.0}, "ECAL Time Window minimum time in ns"};
  Gaudi::Property<float> m_ecalEndcapTimeWindowMax{
      this, "ECALEndcapTimeWindowMax", {100.0}, "ECAL Endcap Time Window maximum time in ns"};
  Gaudi::Property<float> m_ecalBarrelTimeWindowMax{
      this, "ECALBarrelTimeWindowMax", {100.0}, "ECAL Barrel Time Window maximum time in ns"};
  Gaudi::Property<float> m_ecalDeltaTimeHitResolution{
      this, "ECALDeltaTimeHitResolution", {10.0}, "ECAL Minimum Delta Time in ns for resolving two hits"};
  Gaudi::Property<float> m_ecalTimeResolution{
      this, "ECALTimeResolution", {10.0}, "ECAL Time Resolution used to smear hit times"};
  Gaudi::Property<bool> m_ecalSimpleTimingCut{
      this,
      "ECALSimpleTimingCut",
      {true},
      "Use simple time window cut on hit times? If false: use original hit-time clustering algorithm. If true: use "
      "time window defined by ECALBarrelTimeWindowMin and ECALBarrelTimeWindowMax"};

  // timing parameters for HCAL
  Gaudi::Property<int> m_useHcalTiming{this, "UseHcalTiming", {1}, "Use HCAL hit times"};
  Gaudi::Property<int> m_hcalCorrectTimesForPropagation{
      this, "HCALCorrectTimesForPropagation", {0}, "Correct HCAL hit times for propagation: radial distance/c"};
  Gaudi::Property<float> m_hcalTimeWindowMin{this, "HCALTimeWindowMin", {-10.0}, "HCAL Time Window minimum time in ns"};
  Gaudi::Property<float> m_hcalEndcapTimeWindowMax{
      this, "HCALEndcapTimeWindowMax", {100.0}, "HCAL Endcap Time Window maximum time in ns"};
  Gaudi::Property<float> m_hcalBarrelTimeWindowMax{
      this, "HCALBarrelTimeWindowMax", {100.0}, "HCAL Barrel Time Window maximum time in ns"};
  Gaudi::Property<float> m_hcalDeltaTimeHitResolution{
      this, "HCALDeltaTimeHitResolution", {10.0}, "HCAL Minimum Delta Time in ns for resolving two hits"};
  Gaudi::Property<float> m_hcalTimeResolution{
      this, "HCALTimeResolution", {10.0}, "HCAL Time Resolution used to smear hit times"};
  Gaudi::Property<bool> m_hcalSimpleTimingCut{
      this,
      "HCALSimpleTimingCut",
      {true},
      "Use simple time window cut on hit times? If false: use original hit-time clustering algorithm. If true: use "
      "time window defined by HCALBarrelTimeWindowMin and HCALBarrelTimeWindowMax"};

  // parameters for extra ECAL digitization effects
  Gaudi::Property<float> m_calibEcalMip{
      this, "CalibECALMIP", {1.0e-4}, "Calibration to convert ECAL deposited energy to MIPs"};
  Gaudi::Property<int> m_applyEcalDigi{
      this,
      "ECALApplyRealisticDigi",
      {0},
      "Apply realistic digitisation to ECAL hits? (0=none, 1=silicon, 2=scintillator)"};
  Gaudi::Property<float> m_ecal_PPD_pe_per_mip{
      this, "ECAL_PPD_PE_per_MIP", {7.0}, "# photoelectrons per MIP (scintillator): used to poisson smear #PEs if > 0"};
  Gaudi::Property<int> m_ecal_PPD_n_pixels{
      this,
      "ECAL_PPD_N_Pixels",
      {10000},
      "ECAL total number of MPPC/SiPM pixels for implementation of saturation effect"};
  Gaudi::Property<float> m_ecal_misCalibNpix{
      this,
      "ECAL_PPD_N_Pixels_uncertainty",
      {0.05},
      "ECAL fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> m_misCalibEcal_uncorrel{
      this,
      "ECAL_miscalibration_uncorrel",
      {0.0},
      "Uncorrelated ECAL random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<bool> m_misCalibEcal_uncorrel_keep{
      this,
      "ECAL_miscalibration_uncorrel_memorise",
      {false},
      "Store uncorrelated ECAL miscalbrations in memory? (WARNING: Can take a lot of memory if used...)"};
  Gaudi::Property<float> m_misCalibEcal_correl{
      this,
      "ECAL_miscalibration_correl",
      {0.0},
      "Correlated ECAL random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<float> m_deadCellFractionEcal{
      this, "ECAL_deadCellRate", {0.0}, "ECAL random dead cell fraction (as a fraction: [0;1])"};
  Gaudi::Property<bool> m_deadCellEcal_keep{
      this,
      "ECAL_deadCell_memorise",
      {false},
      "Store dead ECAL cells in memory? (WARNING: Can take a lot of memory if used...)"};
  Gaudi::Property<float> m_strip_abs_length{
      this, "ECAL_strip_absorbtionLength", {1000000.0}, "Length scale for absorbtion along scintillator strip (mm)"};
  Gaudi::Property<float> m_ecal_pixSpread{
      this, "ECAL_pixel_spread", {0.05}, "Variation of MPPC/SiPM pixels capacitance in ECAL (as a fraction: 0.01=1%)"};
  Gaudi::Property<float> m_ecal_elec_noise{
      this, "ECAL_elec_noise_mips", {0.0}, "Typical electronics noise for ECAL (in MIP units)"};
  Gaudi::Property<float> m_ehEnergy{
      this, "energyPerEHpair", {3.6}, "Energy required to create e-h pair in silicon (in eV)"};
  Gaudi::Property<float> m_ecalMaxDynMip{
      this, "ECAL_maxDynamicRange_MIP", {2500.0}, "Maximum of electronis dynamic range for ECAL (in MIPs)"};
  Gaudi::Property<int> m_ecalStrip_default_nVirt{
      this, "StripEcal_default_nVirtualCells", {9}, "Default number of virtual cells (used if not found in gear file)"};
  Gaudi::Property<std::string> m_ecal_deafult_layer_config{
      this,
      "ECAL_default_layerConfig",
      {"000000000000000"},
      "Default ECAL layer configuration (used if not found in gear file"};

  // parameters for extra HCAL digitization effects
  Gaudi::Property<float> m_calibHcalMip{
      this, "CalibHCALMIP", {1.0e-4}, "Calibration to convert HCAL deposited energy to MIPs"};
  Gaudi::Property<int> m_applyHcalDigi{
      this,
      "HCALApplyRealisticDigi",
      {0},
      "Apply realistic digitisation to HCAL hits? (0=none, 1=silicon, 2=scintillator)"};
  Gaudi::Property<float> m_hcal_PPD_pe_per_mip{
      this,
      "HCAL_PPD_PE_per_MIP",
      {7.0},
      "# photo-electrons per MIP (scintillator): used to poisson smear #PEs if > 0"};
  Gaudi::Property<int> m_hcal_PPD_n_pixels{
      this,
      "HCAL_PPD_N_Pixels",
      {10000},
      "HCAL total number of MPPC/SiPM pixels for implementation of saturation effect"};
  Gaudi::Property<float> m_hcal_misCalibNpix{
      this,
      "HCAL_PPD_N_Pixels_uncertainty",
      {0.05},
      "HCAL fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> m_misCalibHcal_uncorrel{
      this,
      "HCAL_miscalibration_uncorrel",
      {0.0},
      "Uncorrelated HCAL random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<bool> m_misCalibHcal_uncorrel_keep{
      this,
      "HCAL_miscalibration_uncorrel_memorise",
      {false},
      "Store oncorrelated HCAL miscalbrations in memory? (WARNING: Can take a lot of memory if used...)"};
  Gaudi::Property<float> m_misCalibHcal_correl{
      this,
      "HCAL_miscalibration_correl",
      {0.0},
      "Correlated HCAL random gaussian miscalibration (as a fraction: 1.0 = 100%)"};
  Gaudi::Property<float> m_deadCellFractionHcal{
      this, "HCAL_deadCellRate", {0.0}, "HCAL random dead cell fraction (as a fraction: [0;1])"};
  Gaudi::Property<bool> m_deadCellHcal_keep{
      this,
      "HCAL_deadCell_memorise",
      {false},
      "Store dead HCAL cells in memory? (WARNING: Can take a lot of memory if used...)"};
  Gaudi::Property<float> m_hcal_pixSpread{
      this, "HCAL_pixel_spread", {0.05}, "Variation of MPPC/SiPM pixels capacitance in HCAL (as a fraction: 0.01=1%)"};
  Gaudi::Property<float> m_hcal_elec_noise{
      this, "HCAL_elec_noise_mips", {0.}, "Typical electronics noise for HCAL (in MIP units)"};
  Gaudi::Property<float> m_hcalMaxDynMip{
      this, "HCAL_maxDynamicRange_MIP", {2500.}, "Maximum of electronis dynamic range for HCAL (in MIPs)"};
  //Gaudi::Property<int> m_hcalStrip_default_nVirt{this, "StripHcal_default_nVirtualCells", {9}, "Default number of virtual cells (used if not found in gear file)"};
  //Gaudi::Property<std::string> m_hcal_deafult_layer_config{this, "HCAL_default_layerConfig", {"000000000000000"}, "Default HCAL layer configuration (used if not found in gear file"};
  //Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
  //"The name of the DD4hep constant that contains the Encoding string for tracking detectors"};

  SmartIF<IGeoSvc>         m_geoSvc;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  //const float slop = 0.25; // (mm)

  virtual void fillECALGaps(std::vector<edm4hep::MutableCalorimeterHit*> m_calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS],
                            std::vector<int> m_calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS]) const;

  float digitalEcalCalibCoeff(int layer) const;
  float analogueEcalCalibCoeff(int layer) const;

  float digitalHcalCalibCoeff(CHT::Layout, float energy) const;
  float analogueHcalCalibCoeff(CHT::Layout, int layer) const;

  float ecalEnergyDigi(float energy, int id, float event_correl_miscalib_ecal, CLHEP::MTwistEngine&) const;
  float hcalEnergyDigi(float energy, int id, float event_correl_miscalib_hcal, CLHEP::MTwistEngine&) const;

  float siliconDigi(float energy, CLHEP::MTwistEngine&) const;
  float scintillatorDigi(float energy, bool isEcal, CLHEP::MTwistEngine&) const;

  std::vector<std::pair<int, int>>     getLayerConfig() const;
  void                                 checkConsistency(std::string colName, int layer) const;
  std::pair<int, int>                  getLayerProperties(std::string const& colName, int layer) const;
  int                                  getStripOrientationFromColName(std::string const& colName) const;
  dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag = 0) const;

  float m_zOfEcalEndcap = 0.0;
  float m_barrelPixelSizeT[MAX_LAYERS];
  float m_barrelPixelSizeZ[MAX_LAYERS];
  float m_endcapPixelSizeX[MAX_LAYERS];
  float m_endcapPixelSizeY[MAX_LAYERS];
  float m_barrelStaveDir[MAX_STAVES][2];

  std::unique_ptr<DDScintillatorPpdDigi> m_scEcalDigi{};
  std::unique_ptr<DDScintillatorPpdDigi> m_scHcalDigi{};

  // internal variables
  mutable int m_countWarnings = 0;
  std::string m_ecalLayout    = "";

  std::map<int, float> m_ECAL_cell_miscalibs{};
  std::map<int, bool>  m_ECAL_cell_dead{};
  std::map<int, float> m_HCAL_cell_miscalibs{};
  std::map<int, bool>  m_HCAL_cell_dead{};

  enum { SQUARE, STRIP_ALIGN_ALONG_SLAB, STRIP_ALIGN_ACROSS_SLAB, SIECAL = 0, SCECAL };
};

DECLARE_COMPONENT(DDCaloDigi)
#endif
