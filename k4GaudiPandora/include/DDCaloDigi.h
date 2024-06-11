#ifndef DDCCALODIGI_H
#define DDCCALODIGI_H 1

#include "Gaudi/Property.h"
#include "GaudiAlg/Transformer.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "k4FWCore/BaseClass.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "CalorimeterHitType.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "DDScintillatorPpdDigi.h"
#include "CLHEP/Random/MTwistEngine.h"

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
const int MAX_STAVES =  16;

using retType = std::tuple<
    std::map<std::string, edm4hep::CalorimeterHitCollection>,
    edm4hep::MCRecoCaloAssociationCollection>; //FIXME: Does this need to be a map as well???


struct DDCaloDigi final
  : k4FWCore :: MultiTransformer<
  <retType>(const std::map<std::string, const edm4hep::SimCalorimeterHitCollection&>, 
           const edm4hep::EventHeaderCollection&)> {
  DDCaloDigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::MCRecoCaloAssociationCollection> operator()
  (const edm4hep::SimCalorimeterHitCollection& simCaloHits,
   const edm4hep::EventHeaderCollection& headers) const override;


   private:

  
  Gaudi::Property<float> _thresholdEcal{this, "ECALThreshold", {5.0e-5}, "Threshold for ECAL Hits in GeV"};
  Gaudi::Property<std::string> _unitThresholdEcal{this, "ECALThresholdUnit", {"GeV"}, "Unit for ECAL Threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  Gaudi::Property<float> _thresholdHcal{this, "HCALThreshold", {0.00004}, "Unit for ECAL Threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  Gaudi::Property<std::string> _unitThresholdHcal{this, "HCALThresholdUnit", {"GeV"}, "Unit for HCAL Threshold. Can be \"GeV\", \"MIP\" or \"px\". MIP and px need properly set calibration constants"};
  std::vector<int> ecalLayers;
  ecalLayers.push_back(20);
  ecalLayers.push_back(100);
  Gaudi::Property<std::vector<int>> _ecalLayer{this, "ECALLayers", {ecalLayers}, "Index of ECal Layers"};
  std::vector<int> hcalLayers;
  hcalLayers.push_back(100);
  Gaudi::Property<std::vector<int>> _hcalLayer{this, "ECALLayers", {hcalLayers}, "Index of HCal Layers"};
  std::vector<float> calibrEcal;
  calibrEcal.push_back(40.91);
  calibrEcal.push_back(81.81);
  Gaudi::Property<std::vector<float>> _calibrCoeffEcal{this, "CalibrECAL", {calibrEcal}, "Calibration coefficients for ECAL"};
  std::vector<float> calibrHcalBarrel;
  calibrHcalBarrel.push_back(0.);
  Gaudi::Property<std::vector<float>> _calibrCoeffHcalBarrel{this, "CalibrHCALBarrel", {calibrHcalBarrel}, "Calibration coefficients for Barrel HCAL"};
  std::vector<float> calibrHcalEndCap;
  calibrHcalEndCap.push_back(0.);
  Gaudi::Property<std::vector<float>> _calibrCoeffHcalEndCap{this, "CalibrHCALEndCap", {calibrHcalEndCap}, "Calibration coefficients for EndCap HCAL"};
  std::vector<float> calibrHcalOther;
  calibrHcalOther.push_back(0.);
  Gaudi::Property<std::vector<float>> _calibrCoeffHcalOther{this, "CalibrHCALOther", {calibrHcalOther}, "Calibration coefficients for Other HCAL"};
  Gaudi::Property<int> _digitalEcal{this, "IfDigitalEcal", {0}, "Digital Ecal"};
  Gaudi::Property<int> _mapsEcalCorrection{this, "MapsEcalCorrection", {0}, "Ecal correction for theta dependency of calibration for MAPS"};
  Gaudi::Property<int> _digitalHcal{this, "IfDigitalHcal", {0}, "Digital Hcal"};
  Gaudi::Property<int> _ecalGapCorrection{this, "ECAconst float slop = 0.25; // (mm)LGapCorrection", {1}, "Correct for ECAL gaps"};
  Gaudi::Property<int> _hcalGapCorrection{this, "HCALGapCorrection", {1}, "Correct for HCAL gaps"};
  Gaudi::Property<float> _ecalEndcapCorrectionFactor{this, "ECALEndcapCorrectionFactor", {1.025}, "Energy correction for ECAL endcap"};
  Gaudi::Property<float> _hcalEndcapCorrectionFactor{this, "HCALEndcapCorrectionFactor", {1.025}, "Energy correction for HCAL endcap"};
  Gaudi::Property<float> _ecalGapCorrectionFactor{this, "ECALGapCorrectionFactor", {1.0}, "Factor applied to gap correction"};
  Gaudi::Property<float> _ecalModuleGapCorrectionFactor{this, "ECALModuleGapCorrectionFactor", {0.5}, "Factor applied to module gap correction ECAL"};
  Gaudi::Property<float> _hcalModuleGapCorrectionFactor{this, "HCALModuleGapCorrectionFactor", {0.5}, "Factor applied to module gap correction HCAL"};
  Gaudi::Property<int> _histograms{this, "Histograms", {0}, "Hit times histograms"};
  Gaudi::Property<int> _useEcalTiming{this, "UseEcalTiming", {0}, "Use ECAL hit times"};
  Gaudi::Property<int> _ecalCorrectTimesForPropagation{this, "ECALCorrectTimesForPropagation", {0}, "Correct ECAL hit times for propagation: radial distance/c"};
  Gaudi::Property<float> _ecalTimeWindowMin{this, "ECALTimeWindowMin", {-10.0}, "ECAL Time Window minimum time in ns"};
  Gaudi::Property<float> _ecalEndcapTimeWindowMax{this, "ECALEndcapTimeWindowMax", {100}, "ECAL Endcap Time Window maximum time in ns"};
  Gaudi::Property<float> _ecalBarrelTimeWindowMax{this, "ECALBarrelTimeWindowMax", {100}, "ECAL Barrel Time Window maximum time in ns"};
  Gaudi::Property<float> _ecalDeltaTimeHitResolution{this, "ECALDeltaTimeHitResolution", {10}, "ECAL Minimum Delta Time in ns for resolving two hits"};
  Gaudi::Property<float> _ecalTimeResolution{this, "ECALTimeResolution", {10}, "ECAL Time Resolution used to smear hit times"};
  Gaudi::Property<bool> _ecalSimpleTimingCut{this, "ECALSimpleTimingCut", {true}, "Use simple time window cut on hit times? If false: use original hit-time clustering algorithm. If true: use time window defined by ECALBarrelTimeWindowMin and ECALBarrelTimeWindowMax"};
  Gaudi::Property<int> _useHcalTiming{this, "UseHcalTiming", {1}, "Use HCAL hit times"};
  Gaudi::Property<int> _hcalCorrectTimesForPropagation{this, "HCALCorrectTimesForPropagation", {0}, "Correct HCAL hit times for propagation: radial distance/c"};
  Gaudi::Property<float> _hcalTimeWindowMin{this, "HCALTimeWindowMin", {-10.0}, "HCAL Time Window minimum time in ns"};
  Gaudi::Property<float> _hcalEndcapTimeWindowMax{this, "HCALEndcapTimeWindowMax", {100}, "HCAL Endcap Time Window maximum time in ns"};
  Gaudi::Property<float> _hcalBarrelTimeWindowMax{this, "HCALBarrelTimeWindowMax", {100}, "HCAL Barrel Time Window maximum time in ns"};
  Gaudi::Property<float> _hcalDeltaTimeHitResolution{this, "HCALDeltaTimeHitResolution", {10}, "HCAL Minimum Delta Time in ns for resolving two hits"};
  Gaudi::Property<float> _hcalTimeResolution{this, "HCALTimeResolution", {10}, "HCAL Time Resolution used to smear hit times"};
  Gaudi::Property<bool> _hcalSimpleTimingCut{this, "HCALSimpleTimingCut", {true}, "Use simple time window cut on hit times? If false: use original hit-time clustering algorithm. If true: use time window defined by HCALBarrelTimeWindowMin and HCALBarrelTimeWindowMax"};
  Gaudi::Property<float> _calibEcalMip{this, "CalibECALMIP", {1.0e-4}, "calibration to convert ECAL deposited energy to MIPs"};
  Gaudi::Property<int> _applyEcalDigi{this, "ECAL_apply_realistic_digi", {0}, "apply realistic digitisation to ECAL hits? (0=none, 1=silicon, 2=scintillator)"};
  Gaudi::Property<float> _ecal_PPD_pe_per_mip{this, "ECAL_PPD_PE_per_MIP", {7.0}, "# Photo-electrons per MIP (scintillator): used to poisson smear #PEs if >0"};
  Gaudi::Property<int> _ecal_PPD_n_pixels{this, "ECAL_PPD_N_Pixels",{10000},"ECAL total number of MPPC/SiPM pixels for implementation of saturation effect"};
  Gaudi::Property<float> _ecal_misCalibNpix{this, "ECAL_PPD_N_Pixels_uncertainty", {0.05}, "ECAL fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> _misCalibEcal_uncorrel{this, "ECAL_miscalibration_uncorrel", {0.0}, "uncorrelated ECAL random gaussian miscalibration (as a fraction: 1.0 = 100%"};
  Gaudi::Property<bool> _misCalibEcal_uncorrel_keep{this, "ECAL_miscalibration_uncorrel_memorise", {false}, "store oncorrelated ECAL miscalbrations in memory? (WARNING: can take a lot of memory if used..."};
  Gaudi::Property<float> _misCalibEcal_correl{this, "ECAL_miscalibration_correl", {0.0}, "correlated ECAL random gaussian miscalibration (as a fraction: 1.0 = 100%"};
  Gaudi::Property<float> _deadCellFractionEcal{this, "ECAL_deadCellRate", {0.0}, "ECAL random dead cell fraction (as a fraction: 0->1)"};
  Gaudi::Property<bool> _deadCellEcal_keep{this, "ECAL_deadCell_memorise", {false}, "store dead ECAL cells in memory? (WARNING: can take a lot of memory if used..."};
  Gaudi::Property<float> _strip_abs_length{this, "ECAL_strip_absorbtionLength", {1000000.0}, "length scale for absorbtion along scintillator strip (mm)"};
  Gaudi::Property<float> _ecal_pixSpread{this, "ECAL_pixel_spread", {0.05}, "variation of mppc/sipm pixels capacitance in ECAL (as a fraction: 0.01=1%)"};
  Gaudi::Property<float> _ecal_elec_noise{this, "ECAL_elec_noise_mips", {0.}, "typical electronics noise (ECAL, in MIP units)"};
  Gaudi::Property<float> _ehEnergy{this, "energyPerEHpair", {3.6}, "energy required to create e-h pair in silicon (in eV)"};
  Gaudi::Property<float> _ecalMaxDynMip{this, "ECAL_maxDynamicRange_MIP", {2500.}, "maximum of dynamic range for ECAL (in MIPs)"};
  Gaudi::Property<int> _ecalStrip_default_nVirt{this, "StripEcal_default_nVirtualCells", {9}, "default number of virtual cells (used if not found in gear file)"};
  Gaudi::Property<std::string> _ecal_deafult_layer_config{this, "ECAL_default_layerConfig", {"000000000000000"}, "default ECAL layer configuration (used if not found in gear file"};
  Gaudi::Property<float> _calibHcalMip{this, "CalibHCALMIP", {1.0e-4}, "calibration to convert HCAL deposited energy to MIPs"};
  Gaudi::Property<int> _applyHcalDigi{this, "HCAL_apply_realistic_digi", {0}, "apply realistic digitisation to HCAL hits? (0=none, 1=silicon, 2=scintillator)"};
  Gaudi::Property<float> _hcal_PPD_pe_per_mip{this, "HCAL_PPD_PE_per_MIP", {7.0}, "# Photo-electrons per MIP (scintillator): used to poisson smear #PEs if >0"};
  Gaudi::Property<int> _hcal_PPD_n_pixels{this, "HCAL_PPD_N_Pixels",{10000},"HCAL total number of MPPC/SiPM pixels for implementation of saturation effect"};
  Gaudi::Property<float> _hcal_misCalibNpix{this, "HCAL_PPD_N_Pixels_uncertainty", {0.05}, "HCAL fractional uncertainty of effective total number of MPPC/SiPM pixels"};
  Gaudi::Property<float> _misCalibHcal_uncorrel{this, "HCAL_miscalibration_uncorrel", {0.0}, "uncorrelated HCAL random gaussian miscalibration (as a fraction: 1.0 = 100%"};
  Gaudi::Property<bool> _misCalibHcal_uncorrel_keep{this, "HCAL_miscalibration_uncorrel_memorise", {false}, "store oncorrelated HCAL miscalbrations in memory? (WARNING: can take a lot of memory if used..."};
  Gaudi::Property<float> _misCalibHcal_correl{this, "HCAL_miscalibration_correl", {0.0}, "correlated HCAL random gaussian miscalibration (as a fraction: 1.0 = 100%"};
  Gaudi::Property<float> _deadCellFractionHcal{this, "HCAL_deadCellRate", {0.0}, "HCAL random dead cell fraction (as a fraction: 0->1)"};
  Gaudi::Property<bool> _deadCellHcal_keep{this, "HCAL_deadCell_memorise", {false}, "store dead HCAL cells in memory? (WARNING: can take a lot of memory if used..."};
  Gaudi::Property<float> _hcal_pixSpread{this, "HCAL_pixel_spread", {0.05}, "variation of mppc/sipm pixels capacitance in HCAL (as a fraction: 0.01=1%)"};
  Gaudi::Property<float> _hcal_elec_noise{this, "HCAL_elec_noise_mips", {0.}, "typical electronics noise (ECAL, in MIP units)"};
  Gaudi::Property<float> _hcalMaxDynMip{this, "HCAL_maxDynamicRange_MIP", {2500.}, "maximum of dynamic range for HCAL (in MIPs)"};
  //Gaudi::Property<int> _ecalStrip_default_nVirt{this, "StripEcal_default_nVirtualCells", {9}, "default number of virtual cells (used if not found in gear file)"};
  //Gaudi::Property<std::string> _ecal_deafult_layer_config{this, "ECAL_default_layerConfig", {"000000000000000"}, "default ECAL layer configuration (used if not found in gear file"};
 


  
  const float slop = 0.25; // (mm)
  //DDCaloDigi() ;
 // DDCaloDigi(const DDCaloDigi&) = delete;
 // DDCaloDigi& operator=(const DDCaloDigi&) = delete;

  virtual void fillECALGaps() ;
  
  float digitalHcalCalibCoeff(CHT::Layout,float energy );

  float analogueHcalCalibCoeff(CHT::Layout, int layer );

  float digitalEcalCalibCoeff(int layer );

  float analogueEcalCalibCoeff(int layer );


  float ecalEnergyDigi(float energy,int id);
  float ahcalEnergyDigi(float energy, int id);

  float siliconDigi(float energy);
  float scintillatorDigi(float energy, bool isEcal);
  auto combineVirtualStripCells(auto col, bool isBarrel, int stripOrientation );

  int getNumberOfVirtualCells();
  std::vector < std::pair <int, int> > & getLayerConfig();
  void checkConsistency(std::string colName, int layer);
  std::pair < int, int > getLayerProperties( std::string colName, int layer );
  int getStripOrientationFromColName( std::string colName );


  int _nRun = 0;
  int _nEvt = 0;
  
  //LCFlagImpl _flag{};

  

  std::string _outputRelCollection = "";

  float _thresholdEcal = 5.0e-5;
  std::string _unitThresholdEcal = "GeV";
  std::vector<float> _thresholdHcal{};
  std::string _unitThresholdHcal = "GeV";

  int _digitalEcal = 0;
  int _mapsEcalCorrection = 0;
  int _digitalHcal = 0;

  //bool _ECAL_stripHits;

  std::vector<float> _calibrCoeffEcal{};
  std::vector<float> _calibrCoeffHcalBarrel{};
  std::vector<float> _calibrCoeffHcalEndCap{};
  std::vector<float> _calibrCoeffHcalOther{};

  std::vector<int> _ecalLayers{};
  std::vector<int> _hcalLayers{};

  int _ecalGapCorrection = 1;
  float _ecalGapCorrectionFactor = 1;
  float _ecalModuleGapCorrectionFactor = 0.5;
  float _ecalEndcapCorrectionFactor = 1.025;
  float _hcalEndcapCorrectionFactor = 1.025;
  int   _hcalGapCorrection = 1;
  float _hcalModuleGapCorrectionFactor = 0.5;

  std::vector<edm4hep::MutableCalorimeterHit*> _calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS];
  std::vector<int> _calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS];

  float _zOfEcalEndcap = 0.0;
  float _barrelPixelSizeT[MAX_LAYERS];
  float _barrelPixelSizeZ[MAX_LAYERS];
  float _endcapPixelSizeX[MAX_LAYERS];
  float _endcapPixelSizeY[MAX_LAYERS];
  float _barrelStaveDir[MAX_STAVES][2];
  
  int   _histograms = 0;

  // timing
  int   _useEcalTiming = 0;
  int   _ecalCorrectTimesForPropagation = 0;
  float _ecalTimeWindowMin = -10.0;
  float _ecalBarrelTimeWindowMax = 100.0;
  float _ecalEndcapTimeWindowMax = 100.0;
  float _ecalDeltaTimeHitResolution = 10.0;
  float _ecalTimeResolution = 10.0;
  bool  _ecalSimpleTimingCut = true;

  int   _useHcalTiming = 1;
  int   _hcalCorrectTimesForPropagation = 0;
  float _hcalTimeWindowMin = -10.0;
  float _hcalBarrelTimeWindowMax = 100.0;
  float _hcalEndcapTimeWindowMax = 100.0;
  float _hcalDeltaTimeHitResolution = 10.0;
  float _hcalTimeResolution = 10.0;
  bool  _hcalSimpleTimingCut = true;
  
  std::unique_ptr<DDScintillatorPpdDigi> _scEcalDigi{};
  std::unique_ptr<DDScintillatorPpdDigi> _scHcalDigi{};


  // parameters for extra ECAL digitization effects
  float _calibEcalMip = 1.0e-4;       // MIP calibration factor
  int   _applyEcalDigi = 0;           // which realistic calib to apply
  float _ecal_PPD_pe_per_mip = 7;     // # photoelectrons/MIP for MPPC
  int   _ecal_PPD_n_pixels = 10000;   // # pixels in MPPC
  float _ehEnergy = 3.6;              // energy to create e-h pair in silicon
  float _ecal_misCalibNpix = 0.05;    // miscalibration of # MPPC pixels

  float _misCalibEcal_uncorrel = 0.0; // general ECAL miscalibration (uncorrelated between channels)
  bool  _misCalibEcal_uncorrel_keep = false;// if true, use the same ECAL cell miscalibs in each event (requires more memory)
  float _misCalibEcal_correl = 0.0;     // general ECAL miscalibration (100% uncorrelated between channels)

  float _deadCellFractionEcal = 0.0;  // fraction of random dead channels
  bool  _deadCellEcal_keep = false;   // keep same cells dead between events?

  float _strip_abs_length = 1000000;  // absorption length along strip for non-uniformity modeling
  float _ecal_pixSpread = 0.05;       // relative spread of MPPC pixel signal
  float _ecal_elec_noise = 0;         // electronics noise (as fraction of MIP)
  float _ecalMaxDynMip = 2500;        // electronics dynamic range (in terms of MIPs)
  int _ecalStrip_default_nVirt = 9;   // # virtual cells used in Mokka simulation of strips (if available, this is taken from gear file)
  std::string _ecal_deafult_layer_config ="000000000000000";// ECAL layer configuration (if available, this is taken from gear file)

  // parameters for extra AHCAL digitization effects
  float _calibHcalMip = 1.0e-4;       // MIP calibration factor
  int   _applyHcalDigi = 0;           // which realistic calib to apply
  float _hcal_PPD_pe_per_mip = 10;    // # photoelectrons/MIP for MPPC
  int   _hcal_PPD_n_pixels= 400;      // # pixels in MPPC
  float _hcal_misCalibNpix = 0.05;    // miscalibration of # MPPC pixels

  float _misCalibHcal_uncorrel = 0.0; // general ECAL miscalibration (uncorrelated between channels)
  bool  _misCalibHcal_uncorrel_keep = false; // if true, use the same AHCAL cell miscalibs in each event (requires more memory)
  float _misCalibHcal_correl = 0.0;   // general ECAL miscalibration (100% uncorrelated between channels)

  float _deadCellFractionHcal = 0.0;  // fraction of random dead channels
  bool  _deadCellHcal_keep = false;   // keep same cells dead between events?
  float _hcal_pixSpread = 0.0;        // relative spread of MPPC pixel signal
  float _hcal_elec_noise = 0.0;       // electronics noise (as fraction of MIP)
  float _hcalMaxDynMip = 200;         // electronics dynamic range (in terms of MIPs)



  // internal variables
  std::vector < std::pair <int, int> > _layerTypes {};
  int _strip_virt_cells = 999;
  int _countWarnings = 0;
  std::string _ecalLayout = "";

  float _event_correl_miscalib_ecal = 0.0;
  float _event_correl_miscalib_hcal = 0.0;
  
  CLHEP::MTwistEngine *_randomEngineDeadCellEcal = NULL;
  CLHEP::MTwistEngine *_randomEngineDeadCellHcal = NULL;

  std::map < std::pair <int, int> , float > _ECAL_cell_miscalibs{};
  std::map < std::pair <int, int> , bool > _ECAL_cell_dead{};
  std::map < std::pair <int, int> , float > _HCAL_cell_miscalibs{};
  std::map < std::pair <int, int> , bool > _HCAL_cell_dead{};

  enum {
    SQUARE,
    STRIP_ALIGN_ALONG_SLAB,
    STRIP_ALIGN_ACROSS_SLAB,
    SIECAL=0,
    SCECAL
  };

  
  TH1F* fEcal = NULL;
  TH1F* fHcal = NULL;
  TH1F* fEcalC = NULL;
  TH1F* fHcalC = NULL;
  TH1F* fEcalC1 = NULL;
  TH1F* fHcalC1 = NULL;
  TH1F* fEcalC2 = NULL;
  TH1F* fHcalC2 = NULL;
  TH2F* fHcalCvsE = NULL;
  TH2F* fHcalLayer1 = NULL;
  TH2F* fHcalLayer11 = NULL;
  TH2F* fHcalLayer21 = NULL;
  TH2F* fHcalLayer31 = NULL;
  TH2F* fHcalLayer41 = NULL;
  TH2F* fHcalLayer51 = NULL;
  TH2F* fHcalLayer61 = NULL;
  TH2F* fHcalLayer71 = NULL;
  TH1F* fHcalRLayer1 = NULL;
  TH1F* fHcalRLayer11 = NULL;
  TH1F* fHcalRLayer21 = NULL;
  TH1F* fHcalRLayer31 = NULL;
  TH1F* fHcalRLayer41 = NULL;
  TH1F* fHcalRLayer51 = NULL;
  TH1F* fHcalRLayer61 = NULL;
  TH1F* fHcalRLayer71 = NULL;
  TH1F* fHcalRLayerNorm = NULL;

  TH1F* fEcalRLayerNorm = NULL;
  TH2F* fEcalLayer1 = NULL;
  TH2F* fEcalLayer11 = NULL;
  TH2F* fEcalLayer21 = NULL;
  TH1F* fEcalRLayer1 = NULL;
  TH1F* fEcalRLayer11 = NULL;
  TH1F* fEcalRLayer21 = NULL;

} ;
DECLARE_COMPONENT(DDCaloDigi)
#endif










//  ecalCollections.push_back(std::string("EcalRingCollection"));
  // std::vector<std::string> _ecalCollections; // this is for silicon
  // _ecalCollections.push_back(std::string("EcalBarrelCollection"));
  // _ecalCollections.push_back(std::string("EcalEndcapCollection"));
  //Gaudi::Property<std::vector<std::string>> _ecalCollections{this, "ECALCollections", {_ecalCollections}, "ECAL Collection Names"};
  
  //  std::vector<std::string> _hcalCollections;
  // _hcalCollections.push_back(std::string("HcalBarrelRegCollection"));
  // _hcalCollections.push_back(std::string("HcalEndcapRingsCollection"));
  // _hcalCollections.push_back(std::string("HcalEndcapsCollection"));
 //Gaudi::Property<std::vector<std::string>> _hcalCollections{this, "HCALCollections", {_hcalCollections}, "HCAL Collection Names"};
  
  //  std::vector<std::string> _outputEcalCollections{};
  // _outputEcalCollections.push_back(std::string("ECALBarrel"));
  // _outputEcalCollections.push_back(std::string("ECALEndcap"));
  // _outputEcalCollections.push_back(std::string("ECALOther"));
  
  //  std::vector<std::string> _outputHcalCollections{};
  // _outputHcalCollections.push_back(std::string("HCALBarrel"));
  // _outputHcalCollections.push_back(std::string("HCALEndcap"));
  // _outputHcalCollections.push_back(std::string("HCALOther"));

  // Gaudi::Property<std::string> _outputEcalCollections[0]{this, "output_ECALCollections0", {"ECALBarrel"}, "ECAL Collection of real Hits in Barrel"};
  // Gaudi::Property<std::string> _outputEcalCollections[1]{this, "output_ECALCollections1", {"ECALEndcap"}, "ECAL Collection of real Hits in EndCap"};
  // Gaudi::Property<std::string> _outputEcalCollections[2]{this, "output_ECALCollections2", {"ECALOther"}, "ECAL Collection of real Hits other"};
  
  // Gaudi::Property<std::string> _outputHcalCollections[0]{this, "output_HCALCollections0", {"HCALBarrel"}, "HCAL Collection of real Hits in Barrel"};
  // Gaudi::Property<std::string> _outputHcalCollections[1]{this, "output_HCALCollections1", {"HCALEndcap"}, "HCAL Collection of real Hits in EndCap"};
  // Gaudi::Property<std::string> _outputHcalCollections[2]{this, "output_HCALCollections2", {"HCALOther"}, "HCAL Collection of real Hits other"};
  // Gaudi::Property<std::string> _outputRelCollection{this, "RelationOutputCollection", {"RelationCaloHit"}, "CaloHit Relation Collection"};


