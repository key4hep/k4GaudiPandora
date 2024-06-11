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
#include "DDCaloDigi.h"
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Constants.h"
#include "k4FWCore/BaseClass.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Factories.h"
#include "DDRec/DetectorData.h"
#include "DDRec/DetectorData.h"
#include "DDRec/MaterialManager.h"

using namespace std;
using namespace dd4hep ;
using namespace DDSegmentation ;

// protect against rounding errors
// will not find caps smaller than this
const float slop = 0.25;  // (mm)
//const float pi = acos(-1.0); ///FIXME
//const float twopi = 2.0*pi;  ///FIXME: DD4HEP INTERFERES WITH THESE

DECLARE_COMPONENT(DDCaloDigi)

DDCaloDigi::DDCaloDigi(const std::string& aName, ISvcLocator* aSvcLoc)
    : MultiTransformer(aName, aSvcLoc,
          {
              KeyValues("SimCaloHitCollections", {"SimCalorimeterHitCollection1", "SimCalorimeterHitCollection2",
                                                  "SimCalorimeterHitCollection3"}),
              KeyValues("HeaderName", {"EventHeader"}),
          },
          {KeyValues("output_ECALCollections", {"ECalorimeterHit1", "ECalorimeterHit2", "ECalorimeterHit3"}),
           KeyValues("output_HCALCollections", {"HCalorimeterHit1", "HCalorimeterHit2", "HCalorimeterHit3"}),
           KeyValues("RelationOutputCollection", {"RelationSimCaloHit"})}) {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }
}

//forward declarations. See in DDPandoraPFANewProcessor.cc
dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag = 0);

// // helper struct for string comparision
// struct XToLower{
//   int operator() ( int ch ) {
//     return std::tolower ( ch );
//   }
// }

if (_histograms) {
  fEcal = new TH1F("fEcal", "Ecal time profile", 1000, 0., 1000.0);
  fHcal = new TH1F("fHcal", "Hcal time profile", 1000, 0., 1000.0);

  fEcalC = new TH1F("fEcalC", "Ecal time profile cor", 1000, 0., 1000.0);
  fHcalC = new TH1F("fHcalC", "Hcal time profile cor", 1000, 0., 1000.0);

  fEcalC1 = new TH1F("fEcalC1", "Ecal time profile cor", 100, 0., 1000.0);
  fHcalC1 = new TH1F("fHcalC1", "Hcal time profile cor", 100, 0., 1000.0);

  fEcalC2 = new TH1F("fEcalC2", "Ecal time profile cor", 10, 0., 1000.0);
  fHcalC2 = new TH1F("fHcalC2", "Hcal time profile cor", 10, 0., 1000.0);

  fHcalLayer1     = new TH2F("fHcalLayer1", "Hcal layer 1 map", 300, -4500., 4500.0, 300, -4500, 4500.);
  fHcalLayer11    = new TH2F("fHcalLayer11", "Hcal layer 11 map", 300, -4500., 4500.0, 300, -4500, 4500.);
  fHcalLayer21    = new TH2F("fHcalLayer21", "Hcal layer 21 map", 300, -4500., 4500.0, 300, -4500, 4500.);
  fHcalLayer31    = new TH2F("fHcalLayer31", "Hcal layer 31 map", 300, -4500., 4500.0, 300, -4500, 4500.);
  fHcalLayer41    = new TH2F("fHcalLayer41", "Hcal layer 41 map", 300, -4500., 4500.0, 300, -4500, 4500.);
  fHcalLayer51    = new TH2F("fHcalLayer51", "Hcal layer 51 map", 300, -4500., 4500.0, 300, -4500, 4500.);
  fHcalLayer61    = new TH2F("fHcalLayer61", "Hcal layer 61 map", 300, -4500., 4500.0, 300, -4500, 4500.);
  fHcalLayer71    = new TH2F("fHcalLayer71", "Hcal layer 71 map", 300, -4500., 4500.0, 300, -4500, 4500.);
  fHcalRLayer1    = new TH1F("fHcalRLayer1", "Hcal R layer 1", 50, 0., 5000.0);
  fHcalRLayer11   = new TH1F("fHcalRLayer11", "Hcal R layer 11", 50, 0., 5000.0);
  fHcalRLayer21   = new TH1F("fHcalRLayer21", "Hcal R layer 21", 50, 0., 5000.0);
  fHcalRLayer31   = new TH1F("fHcalRLayer31", "Hcal R layer 31", 50, 0., 5000.0);
  fHcalRLayer41   = new TH1F("fHcalRLayer41", "Hcal R layer 41", 50, 0., 5000.0);
  fHcalRLayer51   = new TH1F("fHcalRLayer51", "Hcal R layer 51", 50, 0., 5000.0);
  fHcalRLayer61   = new TH1F("fHcalRLayer61", "Hcal R layer 61", 50, 0., 5000.0);
  fHcalRLayer71   = new TH1F("fHcalRLayer71", "Hcal R layer 71", 50, 0., 5000.0);
  fHcalRLayerNorm = new TH1F("fHcalRLayerNorm", "Hcal R layer Norm", 50, 0., 5000.0);

  fEcalLayer1     = new TH2F("fEcalLayer1", "Ecal layer 1 map", 1800, -4500., 4500.0, 1800, -4500, 4500.);
  fEcalLayer11    = new TH2F("fEcalLayer11", "Ecal layer 11 map", 1800, -4500., 4500.0, 1800, -4500, 4500.);
  fEcalLayer21    = new TH2F("fEcalLayer21", "Ecal layer 21 map", 1800, -4500., 4500.0, 1800, -4500, 4500.);
  fEcalRLayer1    = new TH1F("fEcalRLayer1", "Ecal R layer 1", 100, 0., 5000.0);
  fEcalRLayer11   = new TH1F("fEcalRLayer11", "Ecal R layer 11", 100, 0., 5000.0);
  fEcalRLayer21   = new TH1F("fEcalRLayer21", "Ecal R layer 21", 100, 0., 5000.0);
  fEcalRLayerNorm = new TH1F("fEcalRLayerNorm", "Ecal R layer Norm", 100, 0., 5000.0);
  m_geoSvc        = serviceLocator()->service(m_geoSvcName);
};

StatusCode DDCaloDigi::initialize() {
  //fHcalCvsE = new TH2F("fHcalCvsE", "Hcal time profile cor",100, 0., 500.0,100,0.,10.)
  _strip_virt_cells = -999;
  _countWarnings    = 0;

  try {
    dd4hep::rec::LayeredCalorimeterData* ecalBarrelData =
        getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

    dd4hep::rec::LayeredCalorimeterData* ecalEndcapData =
        getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& ecalBarrelLayers = ecalBarrelData->layers;
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& ecalEndcapLayers = ecalEndcapData->layers;

    // determine geometry of ECAL
    int symmetry   = ecalBarrelData->inner_symmetry;
    _zOfEcalEndcap = ecalEndcapData->extent[2] / dd4hep::mm;

    // Determine ECAL polygon angles
    // Store radial vectors perpendicular to stave layers in _ecalBarrelStaveDir
    // ASSUMES Mokka Stave numbering 0 = top, then numbering increases anti-clockwise
    if (symmetry > 1) {
      float nFoldSymmetry = static_cast<float>(symmetry);
      float phi0          = ecalBarrelData->phi0 / dd4hep::rad;
      for (int i = 0; i < symmetry; ++i) {
        float phi             = phi0 + i * dd4hep::twopi / nFoldSymmetry;
        _barrelStaveDir[i][0] = cos(phi);
        _barrelStaveDir[i][1] = sin(phi);
      }
    }

    for (unsigned int i = 0; i < ecalBarrelLayers.size(); ++i) {
      _barrelPixelSizeT[i] = ecalBarrelLayers[i].cellSize0;
      _barrelPixelSizeZ[i] = ecalBarrelLayers[i].cellSize1;
      streamlog_out(DEBUG) << "barrel pixel size " << i << " " << _barrelPixelSizeT[i] << " " << _barrelPixelSizeZ[i]
                           << endl;
    }

    for (unsigned int i = 0; i < ecalEndcapLayers.size(); ++i) {
      _endcapPixelSizeX[i] = ecalEndcapLayers[i].cellSize0;
      _endcapPixelSizeY[i] = ecalEndcapLayers[i].cellSize1;
      streamlog_out(DEBUG) << "endcap pixel size " << i << " " << _endcapPixelSizeX[i] << " " << _endcapPixelSizeY[i]
                           << endl;
    }

    _strip_virt_cells = _ecalStrip_default_nVirt;
    streamlog_out(WARNING) << "taking number of virtual cells from steering file (FIXME!): " << _strip_virt_cells
                           << endl;
    _ecalLayout = _ecal_deafult_layer_config;
    streamlog_out(WARNING) << "taking layer layout from steering file (FIXME): " << _ecalLayout << endl;

  } catch (std::exception& e) {
    streamlog_out(ERROR) << "Could not get ECAL parameters from DD4hep!" << endl;
  }

  //convert ECAL thresholds to GeV units
  if (_unitThresholdEcal.compare("GeV") == 0) {
    //ECAL threshold unit is GeV, do nothing
  } else if (_unitThresholdEcal.compare("MIP") == 0) {
    //ECAL threshold unit is MIP, convert via MIP2GeV
    _thresholdEcal *= _calibEcalMip;
  } else if (_unitThresholdEcal.compare("px") == 0) {
    //ECAL threshold unit is pixels, convert via MIP2GeV and lightyield
    _thresholdEcal *= _ecal_PPD_pe_per_mip * _calibEcalMip;
  } else {
    streamlog_out(ERROR) << "could not identify ECAL threshold unit. Please use \"GeV\", \"MIP\" or \"px\"! Aborting."
                         << std::endl;
    assert(0);
  }

  //convert HCAL thresholds to GeV units
  if (_unitThresholdHcal.compare("GeV") == 0) {
    //HCAL threshold unit is GeV, do nothing
  } else if (_unitThresholdHcal.compare("MIP") == 0) {
    //HCAL threshold unit is MIP, convert via MIP2GeV
    for (unsigned int i = 0; i < _thresholdHcal.size(); i++) {
      _thresholdHcal.at(i) *= _calibHcalMip;
    }
  } else if (_unitThresholdHcal.compare("px") == 0) {
    //HCAL threshold unit is pixels, convert via MIP2GeV and lightyield
    for (unsigned int i = 0; i < _thresholdHcal.size(); i++) {
      _thresholdHcal.at(i) *= _hcal_PPD_pe_per_mip * _calibHcalMip;
    }
  } else {
    streamlog_out(ERROR) << "could not identify HCAL threshold unit. Please use \"GeV\", \"MIP\" or \"px\"! Aborting."
                         << std::endl;
    assert(0);
  }

  // set up the scintillator/MPPC digitiser
  _scEcalDigi = std::unique_ptr<DDScintillatorPpdDigi>(new DDScintillatorPpdDigi());
  _scEcalDigi->setPEperMIP(_ecal_PPD_pe_per_mip);
  _scEcalDigi->setCalibMIP(_calibEcalMip);
  _scEcalDigi->setNPix(_ecal_PPD_n_pixels);
  _scEcalDigi->setRandomMisCalibNPix(_ecal_misCalibNpix);
  _scEcalDigi->setPixSpread(_ecal_pixSpread);
  _scEcalDigi->setElecNoise(_ecal_elec_noise);
  _scEcalDigi->setElecRange(_ecalMaxDynMip);
  cout << "ECAL sc digi:" << endl;
  _scEcalDigi->printParameters();

  _scHcalDigi = std::unique_ptr<DDScintillatorPpdDigi>(new DDScintillatorPpdDigi());
  _scHcalDigi->setPEperMIP(_hcal_PPD_pe_per_mip);
  _scHcalDigi->setCalibMIP(_calibHcalMip);
  _scHcalDigi->setNPix(_hcal_PPD_n_pixels);
  _scHcalDigi->setRandomMisCalibNPix(_hcal_misCalibNpix);
  _scHcalDigi->setPixSpread(_hcal_pixSpread);
  _scHcalDigi->setElecNoise(_hcal_elec_noise);
  _scHcalDigi->setElecRange(_hcalMaxDynMip);
  cout << "HCAL sc digi:" << endl;
  _scHcalDigi->printParameters();

  //set up the random engines for ecal and hcal dead cells: (could use a steering parameter though)
  if (_deadCellEcal_keep) {
    _randomEngineDeadCellEcal = new CLHEP::MTwistEngine(0, 0);
  } else {
    _randomEngineDeadCellEcal = 0;
  }

  if (_deadCellHcal_keep) {
    _randomEngineDeadCellHcal = new CLHEP::MTwistEngine(0, 0);
  } else {
    _randomEngineDeadCellHcal = 0;
  }
}

retType DDCaloDigi::operator()(
    const std::map<std::string, const edm4hep::SimCalorimeterHitCollection&>& simCaloHitCollections,
    const edm4hep::EventHeaderCollection&                                     headers) const {
  debug() << " process event : " << headers[0].getEventNumber() << " - run  " << headers[0].getRunNumber()

          << endmsg;  // headers[0].getRunNumber(),headers[0].getEventNumber()

  //auto chschcol = edm4hep::CalorimeterHitCollection(); //collection
  std::map<std::string, edm4hep::CalorimeterHitCollection> outputCollections;
  auto Relcol = edm4hep::MCRecoCaloAssociationCollection();  // Relation collection CalorimeterHit, SimCalorimeterHit

  // copy the flags from the input collection
  //_flag.setBit(LCIO::CHBIT_LONG);
  //_flag.setBit(LCIO::RCHBIT_TIME);  //store timing on output hits.

  // decide on this event's correlated miscalibration
  _event_correl_miscalib_ecal = CLHEP::RandGauss::shoot(1.0, _misCalibEcal_correl);
  _event_correl_miscalib_hcal = CLHEP::RandGauss::shoot(1.0, _misCalibHcal_correl);

  //
  // * Reading Collections of ECAL Simulated Hits *
  //

  for (auto const& inputPair : simCaloHitCollections) {
    std::string colName                = inputPair.first;
    auto&       ecalCol                = outputCollections["new" + colName];
    auto const& inputCaloHitCollection = inputPair.second;  // (or is it ->second ??)
    streamlog_out(DEBUG) << "looking for collection: " << colName << endl;

    if (colName.find("dummy") != string::npos) {
      streamlog_out(DEBUG) << "ignoring input ECAL collection name (looks like dummy name)" << colName << endl;
      continue;
    }

    //fg: need to establish the subdetetcor part here
    //    use collection name as cellID does not seem to have that information
    CHT::Layout caloLayout = layoutFromString(colName);

    try {
      // auto col = evt->getCollection( colName.c_str() ) ;
      const auto initString = cellIDHandle.get();
      //std::string    initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
      dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!
                                                                        // check if decoder contains "layer"
      //CellIDDecoder<SimCalorimeterHit> idDecoder( col );

      // create new collection
      //LCCollectionVec *ecalcol = new LCCollectionVec(LCIO::CALORIMETERHIT);
      //auto ecalcol    = edm4hep::CalorimeterHit();
      //ecalcol->setFlag(_flag.getFlag());

      // if making gap corrections clear the vectors holding pointers to calhits
      if (_ecalGapCorrection != 0) {
        for (int is = 0; is < MAX_STAVES; is++) {
          for (int il = 0; il < MAX_LAYERS; il++) {
            _calHitsByStaveLayer[is][il].clear();
            _calHitsByStaveLayerModule[is][il].clear();
          }
        }
      }

      // deal with strips split into virtual cells
      //  if this collection is a strip which has been split into virtual cells, they need to be recombined
      int orientation = getStripOrientationFromColName(colName);
      if (orientation != SQUARE && getNumberOfVirtualCells() > 1) {
        auto col = combineVirtualStripCells(ecalCol, caloLayout == CHT::barrel, orientation);
      }

      // for (int j(0); j < numElements; ++j) {
      for (const auto& hit : inputCaloHitCollection) {
        //  SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
        float energy = hit.getEnergy();

        // apply threshold cut
        if (energy > _thresholdEcal) {
          int cellid = hit.getCellID();
          int layer  = idDecoder(hit)["layer"];
          int stave  = idDecoder(hit)["stave"];
          int module = idDecoder(hit)["module"];

          // check that layer and assumed layer type are compatible
          checkConsistency(colName, layer);

          // save hits by module/stave/layer if required later
          float calibr_coeff(1.);
          float x    = hit.getPosition()[0];
          float y    = hit.getPosition()[1];
          float z    = hit.getPosition()[2];
          float r    = sqrt(x * x + y * y + z * z);
          float rxy  = sqrt(x * x + y * y);
          float cost = fabs(z) / r;

          if (z > 0 && _histograms) {
            if (layer == 1)
              fEcalLayer1->Fill(x, y);
            if (layer == 11)
              fEcalLayer11->Fill(x, y);
            if (layer == 21)
              fEcalLayer21->Fill(x, y);
            if (layer == 1)
              fEcalRLayer1->Fill(rxy);
            if (layer == 11)
              fEcalRLayer11->Fill(rxy);
            if (layer == 21)
              fEcalRLayer21->Fill(rxy);
          }
          if (_digitalEcal) {
            calibr_coeff = this->digitalEcalCalibCoeff(layer);
            if (_mapsEcalCorrection) {
              if (caloLayout == CHT::barrel) {
                float correction = 1.1387 - 0.068 * cost - 0.191 * cost * cost;
                calibr_coeff /= correction;
              } else {
                float correction = 0.592 + 0.590 * cost;
                calibr_coeff /= correction;
              }
            }
          } else {
            calibr_coeff = this->analogueEcalCalibCoeff(layer);
          }
          // if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff *= _ecalEndcapCorrectionFactor;
          if (caloLayout != CHT::barrel)
            calibr_coeff *= _ecalEndcapCorrectionFactor;  // more robust

          // if you want to understand the timing cut code, please refer to the hcal timing cut further below. it is functionally identical, but has comments, explanations and excuses.
          if (_useEcalTiming) {
            float ecalTimeWindowMax = _ecalEndcapTimeWindowMax;
            if (caloLayout == CHT::barrel)
              ecalTimeWindowMax = _ecalBarrelTimeWindowMax;
            float dt = r / 300. - 0.1;

            std::vector<bool> used(n, false);
            //for(unsigned int i =0; i<n;i++) used[i] = false;
            auto               EsingleHit  = hit.getContributions();
            int                count       = 0;
            float              eCellInTime = 0.;
            float              eCellOutput = 0.;
            const unsigned int n           = EsingleHit.size();
            ;  //number of subhits of this SimHit
            for (unsigned int i_t = 0; i_t < n; i_t++) {
              float timei     = EsingleHit[i_t].getTime();
              float energyi   = EsingleHit[i_t].getEnergy();
              float energySum = 0;

              float deltat = 0;
              if (_ecalCorrectTimesForPropagation)
                deltat = dt;
              if (timei - deltat > _ecalTimeWindowMin && timei - deltat < ecalTimeWindowMax) {
                float ecor = energyi * calibr_coeff;
                eCellInTime += ecor;
              }

              if (!used[i_t]) {
                // merge with other hits?
                used[i_t] = true;
                for (unsigned int j_t = i_t + 1; j_t < n; j_t++) {
                  if (!used[j_t]) {
                    float timej   = EsingleHit[j_t].getTime();
                    float energyj = EsingleHit[j_t].getEnergy();
                    if (_ecalSimpleTimingCut) {
                      float deltat_ij = _ecalCorrectTimesForPropagation ? dt : 0;
                      if (timej - deltat_ij > _ecalTimeWindowMin && timej - deltat_ij < ecalTimeWindowMax) {
                        energySum += energyj;
                        if (timej < timei) {
                          timei = timej;
                        }
                      }
                    } else {
                      float deltat_ij = fabs(timei - timej);
                      if (deltat_ij < _ecalDeltaTimeHitResolution) {
                        if (energyj > energyi)
                          timei = timej;
                        energyi += energyj;
                        used[j_t] = true;
                      }
                    }
                  }
                }
                if (_ecalSimpleTimingCut) {
                  used = vector<bool>(n, true);  // mark everything as used to terminate for loop on next run
                  energyi += energySum;  //fill energySum back into energyi to have rest of loop behave the same.
                }
                if (_digitalEcal) {
                  calibr_coeff = this->digitalEcalCalibCoeff(layer);
                  if (_mapsEcalCorrection) {
                    if (caloLayout == CHT::barrel) {
                      float correction = 1.1387 - 0.068 * cost - 0.191 * cost * cost;
                      calibr_coeff /= correction;
                    } else {
                      float correction = 0.592 + 0.590 * cost;
                      calibr_coeff /= correction;
                    }
                  }
                } else {
                  calibr_coeff = this->analogueEcalCalibCoeff(layer);
                }
                // if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff *= _ecalEndcapCorrectionFactor;
                if (caloLayout != CHT::barrel)
                  calibr_coeff *= _ecalEndcapCorrectionFactor;  // more robust

                if (_histograms) {
                  fEcal->Fill(timei, energyi * calibr_coeff);
                  fEcalC->Fill(timei - dt, energyi * calibr_coeff);
                  fEcalC1->Fill(timei - dt, energyi * calibr_coeff);
                  fEcalC2->Fill(timei - dt, energyi * calibr_coeff);
                }

                // apply extra energy digitisation effects
                energyi = ecalEnergyDigi(energyi, cellID);

                if (energyi > _thresholdEcal) {
                  float timeCor = 0;
                  if (_ecalCorrectTimesForPropagation)
                    timeCor = dt;
                  timei = timei - timeCor;
                  if (timei > _ecalTimeWindowMin && timei < ecalTimeWindowMax) {
                    count++;
                    // CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
                    edm4hep::MutableCalorimeterHit calHit = ecalCol.create();
                    if (_ecalGapCorrection != 0) {
                      _calHitsByStaveLayer[stave][layer].push_back(calhit);
                      _calHitsByStaveLayerModule[stave][layer].push_back(module);
                    }
                    calhit.setCellID(cellid);

                    if (_digitalEcal) {
                      calhit.setEnergy(calibr_coeff);
                    } else {
                      calhit.setEnergy(calibr_coeff * energyi);
                      // calhit->setEnergy(energyi);
                    }

                    eCellOutput += energyi * calibr_coeff;

                    calhit.setTime(timei);
                    calhit.setPosition(hit.getPosition());
                    calhit.setType(CHT(CHT::em, CHT::ecal, caloLayout, layer));
                    calhit.setRawHit(hit);

                    // Set relation with LCRelationNavigator
                    auto caloRel = Relcol.create();
                    caloRel.setRec(calHit);
                    caloRel.setSim(hit);
                    // calohitNav.addRelation(calhit, hit, 1.0);

                  } else {
                    //                    if(caloLayout==CHT::barrel)std::cout << " Drop ECAL Barrel hit : " << timei << " " << calibr_coeff*energyi << std::endl;
                  }
                }
              }
            }
          } else {  // don't use timing
                    // CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
            edm4hep::MutableCalorimeterHit calHit = ecalCol.create();
            if (_ecalGapCorrection != 0) {
              _calHitsByStaveLayer[stave][layer].push_back(calhit);
              _calHitsByStaveLayerModule[stave][layer].push_back(module);
            }
            float energyi = hit.getEnergy();

            // apply extra energy digitisation effects
            energyi = ecalEnergyDigi(energyi, cellid);

            calhit.setCellID(cellid);
            if (_digitalEcal) {
              calhit.setEnergy(calibr_coeff);
            } else {
              calhit.setEnergy(calibr_coeff * energyi);
            }
            calhit.setTime(0);
            calhit.setPosition(hit.getPosition());
            calhit.setType(CHT(CHT::em, CHT::ecal, caloLayout, layer));
            calhit.setRawHit(hit);

            // Set relation with LCRelationNavigator
            auto caloRel = Relcol.create();
            caloRel.setRec(calHit);
            caloRel.setSim(hit);
            //calohitNav.addRelation(calhit, hit, 1.0);
          }  // timing if...else

          //    std::cout << hit->getTimeCont(0) << " count = " << count <<  " E ECAL = " << energyCal << " - " << eCellInTime << " - " << eCellOutput << std::endl;
        }  // energy threshold
      }
      // if requested apply gap corrections in ECAL ?
      if (_ecalGapCorrection != 0)
        this->fillECALGaps();
      // add ECAL collection to event
      // ecalcol->parameters().setValue(LCIO::CellIDEncoding,initString);
      return std::make_tuple(std::move(outputCollections), std::move(caloRel));
      //evt->addCollection(ecalcol,_outputEcalCollections[i].c_str());
    } catch (DataNotAvailableException& e) {
      streamlog_out(DEBUG) << "could not find input ECAL collection " << colName << std::endl;
    }
  }
  //}

  if (_histograms) {
    // fill normalisation of HCAL occupancy plots
    for (float x = 15; x < 3000; x += 30) {
      for (float y = 15; y < 3000; y += 30) {
        if (x > 430 || y > 430) {
          float r = sqrt(x * x + y * y);
          fHcalRLayerNorm->Fill(r, 4.);
        }
      }
    }

    // fill normalisation of ECAL occupancy plots
    for (float x = 2.5; x < 3000; x += 5) {
      for (float y = 2.5; y < 3000; y += 5) {
        float r = sqrt(x * x + y * y);
        if (r > 235)
          fEcalRLayerNorm->Fill(r, 4.);
      }
    }
  }

  //
  // * Reading HCAL Collections of Simulated Hits *
  //

  for (auto const& inputPair : simCaloHitCollections) {
    std::string colName                = inputPair.first;
    auto&       hcalCol                = outputCollections["new" + colName];
    auto const& inputCaloHitCollection = inputPair.second;  // (or is it ->second ??)
    streamlog_out(DEBUG) << "looking for collection: " << colName << endl;

    if (colName.find("dummy") != string::npos) {
      streamlog_out(DEBUG) << "ignoring input HCAL collection name (looks like dummy name)" << colName << endl;
      continue;
    }

    //fg: need to establish the subdetetcor part here
    //    use collection name as cellID does not seem to have that information
    CHT::Layout caloLayout = layoutFromString(colName);

    try {
      //LCCollection * col = evt->getCollection( _hcalCollections[i].c_str() ) ;
      std::string                           initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
      dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!

      hcalcol->setFlag(_flag.getFlag());
      // for (int j(0); j < numElements; ++j) { //loop over all SimCalorimeterHits in this collection
      for (const auto& hit : inputCaloHitCollection) {
        float energy = hit.getEnergy();

        if (energy >
            _thresholdHcal[0] /
                2) {  //preselect for SimHits with energySum>threshold. Doubtful at least, as lower energy hit might fluctuate up and still be counted
          int   cellid = hit.getCellID();
          float calibr_coeff(1.);
          int   layer = idDecoder(hit)["layer"];
          // NOTE : for a digital HCAL this does not allow for varying layer thickness
          // with depth - would need a simple mod to make response proportional to layer thickness
          if (_digitalHcal) {
            calibr_coeff = this->digitalHcalCalibCoeff(caloLayout, energy);
          } else {
            calibr_coeff = this->analogueHcalCalibCoeff(caloLayout, layer);
          }
          // if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff*=_hcalEndcapCorrectionFactor;
          if (caloLayout != CHT::barrel)
            calibr_coeff *= _hcalEndcapCorrectionFactor;  // more robust, is applied to ALL hits outside of barrel.

          //float energyCal = energy*calibr_coeff
          float x = hit.getPosition()[0];
          float y = hit.getPosition()[1];
          float z = hit.getPosition()[2];
          //      float r = sqrt(x*x+y*y);
          if (_useHcalTiming) {
            float hcalTimeWindowMax;
            if (caloLayout == CHT::barrel) {  //current SimHit is in barrel, use barrel timing cut
              hcalTimeWindowMax = _hcalBarrelTimeWindowMax;
            } else {  //current simhit is not in barrel, use endcap timing cut
              hcalTimeWindowMax = _hcalEndcapTimeWindowMax;
            }

            float r = sqrt(
                x * x + y * y +
                z * z);  //this is a crude approximation. assumes initial particle originated at the very center of the detector.
            float              dt         = r / 300 - 0.1;  //magic numbers! ~
            auto               HsingleHit = hit.getContributions();
            const unsigned int n          = HsingleHit.size();
            ;  //number of subhits of this SimHit

            std::vector<bool> used(n, false);

            int count = 0;

            for (unsigned int i_t = 0; i_t < n; i_t++) {      // loop over all subhits
              float timei     = HsingleHit[i_t].getTime();    //absolute hit timing of current subhit
              float energyi   = HsingleHit[i_t].getEnergy();  //energy of current subhit
              float energySum = 0;
              //float deltat = 0;
              //if(_hcalCorrectTimesForPropagation)deltat=dt;
              //deltat now carries hit timing correction.
              //std::cout <<"outer:" << i << " " << n << std::endl;

              //idea of the following section:
              //if simpletimingcut == false
              //sum up hit energies which lie within one calo timing resolution to "timecluster" of current subhit
              //then treat each single timecluster as one hit over threshold and digitise separately. this means there can be more than one CalorimeterHit with the same cellIDs, but different hit times (!)
              //
              //if simpletimingcut == true
              //i'm very sorry. this is the worst code you will ever see.
              //sum up hit energies within timeWindowMin and timeWindowMax, use earliest subhit in this window as hit time for resulting calohit.
              //only one calorimeterhit will be generated from this.

              if (!used
                      [i_t]) {  //if current subhit has not been merged with previous hits already, take current hit as starting point to merge hits
                // merge with other hits?
                used[i_t] = true;
                for (unsigned int j_t = i_t; j_t < n; j_t++) {  //loop through all hits after current hit
                  //std::cout << "inner:" << i << " " << j << " " << n << std::endl;
                  if (!used[j_t]) {
                    float timej   = HsingleHit[j_t].getTime();
                    float energyj = HsingleHit[j_t].getEnergy();
                    //              std::cout << " HCAL  deltat_ij : " << deltat_ij << std::endl;
                    if (_hcalSimpleTimingCut) {
                      float deltat_ij = _hcalCorrectTimesForPropagation ? dt : 0;
                      if (timej - deltat_ij > _hcalTimeWindowMin && timej - deltat_ij < hcalTimeWindowMax) {
                        energySum += energyj;
                        if (timej < timei) {
                          timei = timej;  //use earliest hit time for simpletimingcut
                        }
                      }
                    } else {
                      float deltat_ij = fabs(timei - timej);
                      if (deltat_ij <
                          _hcalDeltaTimeHitResolution) {  //if this subhit is close to current subhit, add this hit's energy to timecluster
                        if (energyj > energyi)
                          timei =
                              timej;  //this is probably not what was intended. i guess this should find the largest hit of one timecluster and use its hittime for the cluster, but instead it compares the current hit energy to the sum of already found hit energies
                        //std::cout << timei << " - " << timej << std::endl;
                        //std::cout << energyi << " - " << energyj << std::endl;
                        energyi += energyj;
                        used[j_t] = true;
                        //std::cout << timei << " " << energyi << std::endl;
                      }
                    }
                  }
                }
                if (_hcalSimpleTimingCut) {
                  used = vector<bool>(
                      n, true);  //mark everything as used to terminate loop. this is worse than goto. i'm sorry.
                  energyi += energySum;  //fill energySum back into energyi to have rest of loop behave the same.
                }

                //variables and their behaviour at this point:
                //if SimpleTimingCut == false
                //energyi carries the sum of subhit energies within +- one hcal time resolution - the timecluster energy.
                //timei carries something vaguely similar to the central hit time of the merged subhits

                //if SimpleTimingCut == true
                //energyi carries the sum of subhit energies within timeWindowMin and timeWindowMax
                //timei carries the time of the earliest hit within this window

                // apply extra energy digitisation effects
                energyi = ahcalEnergyDigi(energyi, cellid);  //this only uses the current subhit "timecluster"!

                if (energyi > _thresholdHcal[0]) {  //now would be the correct time to do threshold comparison
                  float timeCor = 0;
                  if (_hcalCorrectTimesForPropagation)
                    timeCor = dt;
                  timei = timei - timeCor;
                  if (timei > _hcalTimeWindowMin &&
                      timei <
                          hcalTimeWindowMax) {  //if current subhit timecluster is within specified timing window, create new CalorimeterHit and add to collections etc.
                    count++;
                    edm4hep::MutableCalorimeterHit calHit = hcalCol.create();
                    calhit.setCellID(cellid);
                    if (_digitalHcal) {
                      calhit.setEnergy(calibr_coeff);
                      //eCellOutput+= calibr_coeff;
                    } else {
                      calhit.setEnergy(calibr_coeff * energyi);
                      //eCellOutput+= energyi*calibr_coeff;
                    }
                    calhit.setTime(timei);
                    calhit.setPosition(hit.getPosition());
                    calhit.setType(CHT(CHT::had, CHT::hcal, caloLayout, layer));
                    // Set relation with LCRelationNavigator
                    auto hcaloRel = Relcol.create();
                    hcaloRel.setRec(calHit);
                    hcaloRel.setSim(hit);
                    //calohitNav.addRelation(calhit, hit, 1.0);

                  } else {
                    //              std::cout << "Drop HCAL hit : " << timei << " " << calibr_coeff*energyi << std::endl;
                  }
                }
              }
              //std::cout <<"hello" << std::endl;
            }
          } else {  // don't use timing
                    // CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
            edm4hep::MutableCalorimeterHit calHit = hcalCol.create();
            calhit->setCellID(cellid);
            float energyi = hit.getEnergy();

            // apply realistic digitisation
            energyi = ahcalEnergyDigi(energyi, cellid);

            if (_digitalHcal) {
              calhit.setEnergy(calibr_coeff);
            } else {
              calhit.setEnergy(calibr_coeff * energyi);
            }
            //eCellOutput+= energyi*calibr_coeff;
            calhit.setTime(0);
            calhit.setPosition(hit.getPosition());
            calhit.setType(CHT(CHT::had, CHT::hcal, caloLayout, layer));
            calhit.setRawHit(hit);

            auto hcaloRel = Relcol.create();
            hcaloRel.setRec(calHit);
            hcaloRel.setSim(hit);
          }

          // std::cout << hit->getTimeCont(0) << " count = " << count <<  " EHCAL = " << energyCal << " - " << eCellInTime << " - " << eCellOutput << std::endl;
        }
      }
      // add HCAL collection to event
      // hcalcolparameters().setValue(LCIO::CellIDEncoding,initString);
      // evt->addCollection(hcalcol,_outputHcalCollections[i].c_str());
      return std::make_tuple(std::move(outputCollections), std::move(hcaloRel));
    } catch (DataNotAvailableException& e) {
      streamlog_out(DEBUG) << "could not find input HCAL collection " << colName << std::endl;
    }
  }

  // Create and add relation collection for ECAL/HCAL to event
  //chschcol = calohitNav.createLCCollection();
  // evt->addCollection(chschcol,_outputRelCollection.c_str());

  // _nEvt++;
}

StatusCode DDCaloDigi::finalize() { return StatusCode::SUCCESS; }  //fix

if (_histograms) {
  TFile* hfile = new TFile("calTiming.root", "recreate");
  fEcal->TH1F::Write();
  fHcal->TH1F::Write();
  fEcalC->TH1F::Write();
  fHcalC->TH1F::Write();
  fEcalC1->TH1F::Write();
  fHcalC1->TH1F::Write();
  fEcalC2->TH1F::Write();
  fHcalC2->TH1F::Write();
  //fHcalCvsE->TH2F::Write();
  fHcalLayer1->TH2F::Write();
  fHcalLayer11->TH2F::Write();
  fHcalLayer21->TH2F::Write();
  fHcalLayer31->TH2F::Write();
  fHcalLayer41->TH2F::Write();
  fHcalLayer51->TH2F::Write();
  fHcalLayer61->TH2F::Write();
  fHcalLayer71->TH2F::Write();
  fHcalRLayer1->TH1F::Write();
  fHcalRLayer11->TH1F::Write();
  fHcalRLayer21->TH1F::Write();
  fHcalRLayer31->TH1F::Write();
  fHcalRLayer41->TH1F::Write();
  fHcalRLayer51->TH1F::Write();
  fHcalRLayer61->TH1F::Write();
  fHcalRLayer71->TH1F::Write();
  fHcalRLayerNorm->TH1F::Write();
  fEcalLayer1->TH2F::Write();
  fEcalLayer11->TH2F::Write();
  fEcalLayer21->TH2F::Write();
  fEcalRLayer1->TH1F::Write();
  fEcalRLayer11->TH1F::Write();
  fEcalRLayer21->TH1F::Write();
  fEcalRLayerNorm->TH1F::Write();

  hfile->Close();
  delete hfile;
}

//delete randomengines if needed
if (_randomEngineDeadCellHcal != 0) {
  delete _randomEngineDeadCellHcal;
}

if (_randomEngineDeadCellEcal != 0) {
  delete _randomEngineDeadCellEcal;
}

void DDCaloDigi::fillECALGaps() {
  // Loop over hits in the Barrel
  // For each layer calculated differences in hit positions
  // Look for gaps based on expected separation of adjacent hits
  // loop over staves and layers

  for (int is = 0; is < MAX_STAVES; ++is) {
    for (int il = 0; il < MAX_LAYERS; ++il) {
      if (_calHitsByStaveLayer[is][il].size() > 1) {
        // compare all pairs of hits just once (j>i)

        for (unsigned int i = 0; i < _calHitsByStaveLayer[is][il].size() - 1; ++i) {
          edm4hep::MutableCalorimeterHit hiti    = _calHitsByStaveLayer[is][il][i];
          int                            modulei = _calHitsByStaveLayerModule[is][il][i];
          float                          xi      = hiti.getPosition()[0];
          float                          yi      = hiti.getPosition()[1];
          float                          zi      = hiti.getPosition()[2];

          for (unsigned int j = i + 1; j < _calHitsByStaveLayer[is][il].size(); ++j) {
            edm4hep::MutableCalorimeterHit hitj    = _calHitsByStaveLayer[is][il][j];
            int                            modulej = _calHitsByStaveLayerModule[is][il][j];
            float                          xj      = hitj.getPosition()[0];
            float                          yj      = hitj.getPosition()[1];
            float                          zj      = hitj.getPosition()[2];
            float                          dz      = fabs(zi - zj);
            // *** BARREL CORRECTION ***
            if (fabs(zi) < _zOfEcalEndcap && fabs(zj) < _zOfEcalEndcap) {
              // account for stave directions using normals
              // calculate difference in hit postions in z and along stave
              float dx = xi - xj;
              float dy = yi - yj;
              //float dt = fabs(dx*_barrelStaveDir[is][0] + dy*_barrelStaveDir[is][1]);
              float dt = sqrt(dx * dx + dy * dy);  // more generic
              // flags for evidence for gaps
              bool zgap  = false;  // in z direction
              bool tgap  = false;  // along stave
              bool ztgap = false;  // in both z and along stave
              bool mgap  = false;  // gaps between ECAL modules

              // criteria gaps in the z and t direction
              float zminm = 1.0 * _barrelPixelSizeZ[il] - slop;
              float zmin  = 1.0 * _barrelPixelSizeZ[il] + slop;
              float zmax  = 2.0 * _barrelPixelSizeZ[il] - slop;
              float tminm = 1.0 * _barrelPixelSizeT[il] - slop;
              float tmin  = 1.0 * _barrelPixelSizeT[il] + slop;
              float tmax  = 2.0 * _barrelPixelSizeT[il] - slop;

              // criteria for gaps
              // WOULD BE BETTER TO USE GEAR TO CHECK GAPS ARE OF EXPECTED SIZE
              if (dz > zmin && dz < zmax && dt < tminm)
                zgap = true;
              if (dz < zminm && dt > tmin && dt < tmax)
                tgap = true;
              if (dz > zmin && dz < zmax && dt > tmin && dt < tmax)
                ztgap = true;

              if (modulei != modulej) {
                if (dz > zmin && dz < 3.0 * _barrelPixelSizeZ[il] - slop && dt < tmin)
                  mgap = true;
              }

              // found a gap now apply a correction based on area of gap/area of pixel
              if (zgap || tgap || ztgap || mgap) {
                float ecor = 1.;
                float f    = _ecalGapCorrectionFactor;  // fudge
                if (mgap)
                  f = _ecalModuleGapCorrectionFactor;
                if (zgap || mgap)
                  ecor = 1. + f * (dz - _barrelPixelSizeZ[il]) / 2. / _barrelPixelSizeZ[il];
                if (tgap)
                  ecor = 1. + f * (dt - _barrelPixelSizeT[il]) / 2. / _barrelPixelSizeT[il];
                if (ztgap)
                  ecor = 1. + f * (dt - _barrelPixelSizeT[il]) * (dz - _barrelPixelSizeZ[il]) / 4. /
                                  _barrelPixelSizeT[il] / _barrelPixelSizeZ[il];
                float ei = hiti.getEnergy() * ecor;
                float ej = hitj.getEnergy() * ecor;
                hiti.setEnergy(ei);
                hitj.setEnergy(ej);
              }

              // *** ENDCAP CORRECTION ***
            } else if (fabs(zi) > _zOfEcalEndcap && fabs(zj) > _zOfEcalEndcap && dz < 100) {
              float dx    = fabs(xi - xj);
              float dy    = fabs(yi - yj);
              bool  xgap  = false;
              bool  ygap  = false;
              bool  xygap = false;
              // criteria gaps in the z and t direction

              // x and y need to be swapped in different staves of endcap.
              float pixsizex, pixsizey;
              if (is % 2 == 1) {
                pixsizex = _endcapPixelSizeY[il];
                pixsizey = _endcapPixelSizeX[il];
              } else {
                pixsizex = _endcapPixelSizeX[il];
                pixsizey = _endcapPixelSizeY[il];
              }

              float xmin  = 1.0 * pixsizex + slop;
              float xminm = 1.0 * pixsizex - slop;
              float xmax  = 2.0 * pixsizex - slop;
              float ymin  = 1.0 * pixsizey + slop;
              float yminm = 1.0 * pixsizey - slop;
              float ymax  = 2.0 * pixsizey - slop;
              // look for gaps
              if (dx > xmin && dx < xmax && dy < yminm)
                xgap = true;
              if (dx < xminm && dy > ymin && dy < ymax)
                ygap = true;
              if (dx > xmin && dx < xmax && dy > ymin && dy < ymax)
                xygap = true;

              if (xgap || ygap || xygap) {
                // cout <<"NewLDCCaloDigi found endcap gap, adjusting energy! " << xgap << " " << ygap << " " << xygap << " , " << il << endl;
                // cout << "stave " << is <<  " layer " << il << endl;
                // cout << "  dx, dy " << dx<< " " << dy << " , sizes = " << pixsizex << " " << pixsizey << endl;
                // cout << " xmin... " << xmin << " " << xminm << " " << xmax << " ymin... " << ymin << " " << yminm << " " << ymax << endl;

                // found a gap make correction
                float ecor = 1.;
                float f    = _ecalGapCorrectionFactor;  // fudge
                if (xgap)
                  ecor = 1. + f * (dx - pixsizex) / 2. / pixsizex;
                if (ygap)
                  ecor = 1. + f * (dy - pixsizey) / 2. / pixsizey;
                if (xygap)
                  ecor = 1. + f * (dx - pixsizex) * (dy - pixsizey) / 4. / pixsizex / pixsizey;

                // cout << "correction factor = " << ecor << endl;

                hiti.setEnergy(hiti.getEnergy() * ecor);
                hitj.setEnergy(hitj.getEnergy() * ecor);
              }
            }
          }
        }
      }
    }
  }

  return;
}

float DDCaloDigi::digitalHcalCalibCoeff(CHT::Layout caloLayout, float energy) {
  float        calib_coeff = 0;
  unsigned int ilevel      = 0;
  for (unsigned int ithresh = 1; ithresh < _thresholdHcal.size(); ithresh++) {
    // Assume!!!  hit energies are stored as floats, i.e. 1, 2 or 3
    if (energy > _thresholdHcal[ithresh])
      ilevel = ithresh;  // ilevel = 0 , 1, 2
  }

  switch (caloLayout) {
    case CHT::barrel:
      if (ilevel > _calibrCoeffHcalBarrel.size() - 1) {
        streamlog_out(ERROR) << " Semi-digital level " << ilevel
                             << " greater than number of HCAL Calibration Constants (" << _calibrCoeffHcalBarrel.size()
                             << ")" << std::endl;
      } else {
        calib_coeff = _calibrCoeffHcalBarrel[ilevel];
      }
      break;
    case CHT::endcap:
      if (ilevel > _calibrCoeffHcalEndCap.size() - 1) {
        streamlog_out(ERROR) << " Semi-digital level " << ilevel
                             << " greater than number of HCAL Calibration Constants (" << _calibrCoeffHcalEndCap.size()
                             << ")" << std::endl;
      } else {
        calib_coeff = _calibrCoeffHcalEndCap[ilevel];
      }
      break;
    case CHT::plug:
      if (ilevel > _calibrCoeffHcalOther.size() - 1) {
        streamlog_out(ERROR) << " Semi-digital level " << ilevel
                             << " greater than number of HCAL Calibration Constants (" << _calibrCoeffHcalOther.size()
                             << ")" << std::endl;
      } else {
        calib_coeff = _calibrCoeffHcalOther[ilevel];
      }
      break;
    default:
      streamlog_out(ERROR) << " Unknown HCAL Hit Type " << std::endl;
      break;
  }

  return calib_coeff;
}

float DDCaloDigi::analogueHcalCalibCoeff(CHT::Layout caloLayout, int layer) {
  float calib_coeff = 0;

  for (unsigned int k(0); k < _hcalLayers.size(); ++k) {
    int min, max;
    if (k == 0)
      min = 0;
    else
      min = _hcalLayers[k - 1];

    max = _hcalLayers[k];
    if (layer >= min && layer < max) {
      switch (caloLayout) {
        case CHT::barrel:
          calib_coeff = _calibrCoeffHcalBarrel[k];
          break;
        case CHT::endcap:
          calib_coeff = _calibrCoeffHcalEndCap[k];
          break;
        case CHT::plug:
        case CHT::ring:
          calib_coeff = _calibrCoeffHcalOther[k];
          break;
        default:
          streamlog_out(ERROR) << " Unknown HCAL Hit Type " << std::endl;
          break;
      }
    }
  }

  return calib_coeff;
}

float DDCaloDigi::digitalEcalCalibCoeff(int layer) {
  float calib_coeff = 0;

  for (unsigned int k(0); k < _ecalLayers.size(); ++k) {
    int min, max;
    if (k == 0)
      min = 0;
    else
      min = _ecalLayers[k - 1];

    max = _ecalLayers[k];
    if (layer >= min && layer < max) {
      calib_coeff = _calibrCoeffEcal[k];
      break;
    }
  }

  return calib_coeff;
}

float DDCaloDigi::analogueEcalCalibCoeff(int layer) {
  float calib_coeff = 0;

  // retrieve calibration constants
  for (unsigned int k(0); k < _ecalLayers.size(); ++k) {
    int min, max;
    if (k == 0) {
      min = 0;
    } else {
      min = _ecalLayers[k - 1];
    }
    max = _ecalLayers[k];
    if (layer >= min && layer < max) {
      calib_coeff = _calibrCoeffEcal[k];
      break;
    }
  }

  return calib_coeff;
}

float DDCaloDigi::ecalEnergyDigi(float energy, int id) {
  // some extra digi effects (daniel)
  // controlled by _applyEcalDigi = 0 (none), 1 (silicon), 2 (scintillator)

  // small update for time-constant uncorrelated miscalibrations. DJ, Jan 2015

  float e_out(energy);
  if (_applyEcalDigi == 1)
    e_out = siliconDigi(energy);  // silicon digi
  else if (_applyEcalDigi == 2)
    e_out = scintillatorDigi(energy, true);  // scintillator digi

  // add electronics dynamic range
  // Sept 2015: Daniel moved this to the ScintillatorDigi part, so it is applied before unfolding of sipm response
  // if (_ecalMaxDynMip>0) e_out = min (e_out, _ecalMaxDynMip*_calibEcalMip);

  // random miscalib
  if (_misCalibEcal_uncorrel > 0) {
    float miscal(0);
    if (_misCalibEcal_uncorrel_keep) {
      int id;
      if (_ECAL_cell_miscalibs.find(id) !=
          _ECAL_cell_miscalibs.end()) {  // this cell was previously seen, and a miscalib stored
        miscal = _ECAL_cell_miscalibs[id];
      } else {  // we haven't seen this one yet, get a miscalib for it
        miscal                   = CLHEP::RandGauss::shoot(1.0, _misCalibEcal_uncorrel);
        _ECAL_cell_miscalibs[id] = miscal;
      }
    } else {
      miscal = CLHEP::RandGauss::shoot(1.0, _misCalibEcal_uncorrel);
    }
    e_out *= miscal;
  }

  if (_misCalibEcal_correl > 0)
    e_out *= _event_correl_miscalib_ecal;

  // random cell kill
  if (_deadCellFractionEcal > 0) {
    if (_deadCellEcal_keep == true) {
      int id;

      if (_ECAL_cell_dead.find(id) != _ECAL_cell_dead.end()) {  // this cell was previously seen
        if (_ECAL_cell_dead[id] == true) {
          e_out = 0;
        }
      } else {  // we haven't seen this one yet, get a miscalib for it
        bool thisDead       = (CLHEP::RandFlat::shoot(_randomEngineDeadCellEcal, .0, 1.0) < _deadCellFractionEcal);
        _ECAL_cell_dead[id] = thisDead;
        if (thisDead == true) {
          e_out = 0;
        }
      }

    } else {
      if (CLHEP::RandFlat::shoot(0.0, 1.0) < _deadCellFractionEcal)
        e_out = 0;
    }
  }

  return e_out;
}

float DDCaloDigi::ahcalEnergyDigi(float energy, int id) {
  // some extra digi effects (daniel)
  // controlled by _applyHcalDigi = 0 (none), 1 (scintillator/SiPM)

  // small update for time-constant uncorrelated miscalibrations. DJ, Jan 2015

  float e_out(energy);
  if (_applyHcalDigi == 1)
    e_out = scintillatorDigi(energy, false);  // scintillator digi

  // add electronics dynamic range
  // Sept 2015: Daniel moved this to the ScintillatorDigi part, so it is applied before unfolding of sipm response
  //  if (_hcalMaxDynMip>0) e_out = min (e_out, _hcalMaxDynMip*_calibHcalMip);

  // random miscalib
  //  if (_misCalibHcal_uncorrel>0) e_out*=CLHEP::RandGauss::shoot( 1.0, _misCalibHcal_uncorrel );
  if (_misCalibHcal_uncorrel > 0) {
    float miscal(0);
    if (_misCalibHcal_uncorrel_keep) {
      int id;
      if (_HCAL_cell_miscalibs.find(id) !=
          _HCAL_cell_miscalibs.end()) {  // this cell was previously seen, and a miscalib stored
        miscal = _HCAL_cell_miscalibs[id];
      } else {  // we haven't seen this one yet, get a miscalib for it
        miscal                   = CLHEP::RandGauss::shoot(1.0, _misCalibHcal_uncorrel);
        _HCAL_cell_miscalibs[id] = miscal;
      }
    } else {
      miscal = CLHEP::RandGauss::shoot(1.0, _misCalibHcal_uncorrel);
    }
    e_out *= miscal;
  }

  if (_misCalibHcal_correl > 0)
    e_out *= _event_correl_miscalib_hcal;

  // random cell kill
  if (_deadCellFractionHcal > 0) {
    if (_deadCellHcal_keep == true) {
      int id;

      if (_HCAL_cell_dead.find(id) != _HCAL_cell_dead.end()) {  // this cell was previously seen
        if (_HCAL_cell_dead[id] == true) {
          e_out = 0;
        }
      } else {  // we haven't seen this one yet, get a miscalib for it
        bool thisDead       = (CLHEP::RandFlat::shoot(_randomEngineDeadCellHcal, .0, 1.0) < _deadCellFractionHcal);
        _HCAL_cell_dead[id] = thisDead;
        if (thisDead == true) {
          e_out = 0;
        }
      }

    } else {
      if (CLHEP::RandFlat::shoot(0.0, 1.0) < _deadCellFractionHcal)
        e_out = 0;
    }
  }
  return e_out;
}

float DDCaloDigi::siliconDigi(float energy) {
  // applies extra digitisation to silicon hits

  // calculate #e-h pairs
  float nehpairs = 1e9 * energy / _ehEnergy;  // check units of energy! _ehEnergy is in eV, energy in GeV

  // fluctuate it by Poisson
  float smeared_energy = energy * CLHEP::RandPoisson::shoot(nehpairs) / nehpairs;

  // limited electronics dynamic range // Daniel moved electronics dyn range to here
  if (_ecalMaxDynMip > 0)
    smeared_energy = std::min(smeared_energy, _ecalMaxDynMip * _calibEcalMip);

  // add electronics noise
  if (_ecal_elec_noise > 0)
    smeared_energy += CLHEP::RandGauss::shoot(0, _ecal_elec_noise * _calibEcalMip);

  return smeared_energy;
}

float DDCaloDigi::scintillatorDigi(float energy, bool isEcal) {
  // this applies some extra digitisation to scintillator+PPD hits (PPD=SiPM, MPPC)
  // - poisson fluctuates the number of photo-electrons according to #PEs/MIP
  // - applies PPD saturation according to #pixels
  // Daniel Jeans, Jan/Feb 2014.

  float digiEn(0);
  if (isEcal) {
    digiEn = _scEcalDigi->getDigitisedEnergy(energy);  //CHECK!!!
  } else {
    digiEn = _scHcalDigi->getDigitisedEnergy(energy);
  }
  return digiEn;
}

auto DDCaloDigi::combineVirtualStripCells(auto col, bool isBarrel, int stripOrientation) {
  // combines the virtual cells in a strip
  // input collection is for virtual cells
  // returns collection of strips

  // daniel jeans, jan/feb 2014

  // sanity check
  if (stripOrientation == SQUARE) {
    streamlog_out(ERROR) << "DDCaloDigi::combineVirtualStripCells trying to deal with silicon strips??? I do not know "
                            "how to do that, refusing to do anything!"
                         << std::endl;
    return NULL;
  }

  //CellIDDecoder<SimCalorimeterHit> idDecoder( col );
  // This variable is for tracking detectors only, not for calorimeters
  // std::string initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());

  // get the encoding string from the original collection
  const auto                            initString = cellIDHandle.get();
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);

  // make the new output collection
  //edm4hep::MutableSimCalorimeterHit tempstripcol = inputCaloHitCollection.create();;
  auto stripcol = edm4hep::SimCalorimeterHitCollection();
  stripcol.setFlag(_flag.getFlag());
  //FIXME: add initString/ cellIDEncoding to stripcol

  // link between strip CellID and hits storage for strip hits
  std::map<int, edm4hep::MutableSimCalorimeterHit*> newhits;  //SORT

  // loop over input collection
  int numElements = col.size();

  // // sum energy for check
  float tempenergysum(0);
  for (const auto& hit : inputCaloHitCollection) {
    tempenergysum += hit.getEnergy();
  }

  float scTVirtLengthBar(-99);
  float scLVirtLengthBar(-99);
  float scTVirtLengthEnd(-99);
  float scLVirtLengthEnd(-99);

  for (const auto& hit : inputCaloHitCollection) {  // loop over virtual cells

    // for now, ignore preshower layer.
    //   in "usual" (non-preshower) ecal collections, km1 = 0 is the first layer AFTER the preshower
    int km1 = bitFieldCoder(hit)["layer"];  // FIXME

    int gearlayer = km1 + 1;  // in "gearlayer", preshower is 0, first layer after absorber is 1

    // after fix of MOKKA, this should be not needed
    // int gearlayer(-99); // this is the layer index which the GEAR geometry knows about. depends on driver version...
    // if ( _ecal_driver_version == 0 ) { // SEcal04 driver
    // } else if ( _ecal_driver_version == 1 ) { // SEcal05 driver
    //   gearlayer = km1+1;
    // } else {
    //   streamlog_out ( ERROR ) << "ERROR unknown value for ecal driver version, _ecal_driver_version=" << _ecal_driver_version << endl;
    // }
    // assert (gearlayer>=0);

    int m       = bitFieldCoder(hit)["module"];
    int s       = bitFieldCoder(hit)["stave"];
    int i_index = bitFieldCoder(hit)["x"];
    int j_index = bitFieldCoder(hit)["y"];

    //cout << "ORIG: " << km1 << " " << m << " " << s << " " << i_index << " " << j_index << " : " <<
    //  hit->getPosition()[0] << " "<< hit->getPosition()[1] << " "<< hit->getPosition()[2] << endl;

    // should we combine virtual cells in the I or J direction?
    //
    // in barrel:
    //   i is along slab (phi in barrel)
    //   j is across slab (z in barrel)
    //   increasing j = increasing z; increasing i = decreasing phi
    //
    // in endcap:
    //     i is across slab
    //     j is along slab
    //     different staves are rotated around beam axis.
    //        in +ve z endcap (m==6),
    //          stave 0 is at -ve x and +ve y. in this stave, i is parallel to -x, j to +y
    //          stave 3 is at +ve x and +ve y. in this stave, i is parallel to +y, j to +x
    //        -ve z endcap (m==0) is same, just everything flipped in x.

    bool collateAlongI =
        isBarrel ? stripOrientation == STRIP_ALIGN_ALONG_SLAB :  // I is the index along the slab (R-phi in barrel)
            stripOrientation == STRIP_ALIGN_ACROSS_SLAB;         // for some reason seems to be reversed in endcap...!!

    streamlog_out(DEBUG) << "isBarrel = " << isBarrel << " : stripOrientation= " << stripOrientation
                         << " gear layer = " << gearlayer << endl;

    // let's get the length of the virtual cell in the direction of the strip
    //    fixed rare crashes, and streamlined code to get virtualCellSizeAlongStrip. Sept 2014
    float virtualCellSizeAlongStrip(0);

    float* layerpixsize(0);
    float* savedPixSize(0);
    if (isBarrel) {
      if (stripOrientation == STRIP_ALIGN_ACROSS_SLAB) {
        layerpixsize = _barrelPixelSizeT;
        savedPixSize = &scTVirtLengthBar;
      } else {
        layerpixsize = _barrelPixelSizeZ;
        savedPixSize = &scLVirtLengthBar;
      }
    } else {
      if (stripOrientation == STRIP_ALIGN_ACROSS_SLAB) {
        layerpixsize = _endcapPixelSizeX;
        savedPixSize = &scTVirtLengthEnd;
      } else {
        layerpixsize = _endcapPixelSizeY;
        savedPixSize = &scLVirtLengthEnd;
      }
    }

    virtualCellSizeAlongStrip = layerpixsize[gearlayer];
    if (virtualCellSizeAlongStrip > 0) {  // looks OK
      if (*savedPixSize < 0) {            // no saved size yet, keep this one
        *savedPixSize = virtualCellSizeAlongStrip;
      }
    } else {                    // no valid info in gear file for this layer
      if (*savedPixSize > 0) {  // use saved value from other layer
        virtualCellSizeAlongStrip = *savedPixSize;
      } else {  // look at previous layers for one with same orientation (extra check, sept 2014)
        std::vector<std::pair<int, int>> layerTypes = getLayerConfig();

        streamlog_out(DEBUG) << "could not get valid info from gear file..." << endl
                             << "looking through previous layers to get a matching orientation" << endl
                             << "this gear layer " << gearlayer << " type: " << layerTypes[gearlayer].first << " "
                             << layerTypes[gearlayer].second << endl;

        for (int il = gearlayer - 1; il >= 0; il--) {
          // layer types include the preshower at posn "0"
          // gearlayer has preshower as layer 0
          if (layerTypes[il] == layerTypes[gearlayer]) {  // found layer with same setup
            streamlog_out(DEBUG) << "found a match! " << il << " " << layerTypes[il].first << " "
                                 << layerTypes[il].second << " : " << layerpixsize[il] << endl;
            virtualCellSizeAlongStrip = layerpixsize[il];
            *savedPixSize             = virtualCellSizeAlongStrip;
            break;
          }
        }
      }
    }
    streamlog_out(DEBUG) << "virtualCellSizeAlongStrip = " << virtualCellSizeAlongStrip << endl;

    // calculate the new strip's I,J coordinates
    int i_new = collateAlongI ? i_index / getNumberOfVirtualCells() : i_index;
    int j_new = collateAlongI ? j_index : j_index / getNumberOfVirtualCells();

    // apply position-dependent energy correction
    //
    // relative virtual cell position within strip (0 = centre)
    //    distance between centre of virtual cell and centre of strip
    int   relativeIndx = collateAlongI ? i_index : j_index;
    float relativePos =
        relativeIndx % getNumberOfVirtualCells();  // this is index within strip wrt strip end (0->_strip_virt_cells-1)
    relativePos -= (getNumberOfVirtualCells() - 1) / 2.0;  // index wrt strip centre
    relativePos *= virtualCellSizeAlongStrip;              // distance wrt strip centre

    // effect of absorbtion length within scintillator
    //     TODO: should check the polarity is consistent with mppc position, to make sure larger response nearer to mppc....
    float energy_new = hit.getEnergy();

    float energyNonuniformityScaling(1.);
    if (_strip_abs_length > 0) {
      energyNonuniformityScaling = exp(-relativePos / _strip_abs_length);
      energy_new *= energyNonuniformityScaling;
    }

    // calculate ID of the new Strip // FIXME???
    bitFieldCoder["system"] = m;  // FIXME ????
    bitFieldCoder["module"] = m;
    bitFieldCoder["stave"]  = s;
    bitFieldCoder["layer"]  = km1;
    bitFieldCoder["x"]      = i_new;
    bitFieldCoder["y"]      = j_new;
    int cellID              = bitFieldCoder.getCellID;  // FIXME

    // calculate position of centre of this new strip
    // it depends if we are in the barrel or endcap, and in which stave
    //   (n.b. we could do this later, after checking for duplicates. however, keeping it here allows us to
    //            check that we haven't supplied wrong nVirtualCells parameter)

    float stripPos[3] = {0};
    if (isBarrel) {
      if (stripOrientation == STRIP_ALIGN_ACROSS_SLAB) {  // it's along z
        stripPos[0] = hit.getPosition()[0];
        stripPos[1] = hit.getPosition()[1];
        stripPos[2] = hit.getPosition()[2] - relativePos;  // increasing j = increasing z
      } else {                                             // it's along t
        float phi   = atan2(_barrelStaveDir[s][1], _barrelStaveDir[s][0]);
        stripPos[0] = hit.getPosition()[0] - relativePos * cos(phi);  // negative: increasing I = decreasing phi
        stripPos[1] = hit.getPosition()[1] - relativePos * sin(phi);
        stripPos[2] = hit.getPosition()[2];
      }
    } else {  // endcap
      stripPos[0] = hit.getPosition()[0];
      stripPos[1] = hit.getPosition()[1];
      stripPos[2] = hit.getPosition()[2];  // z never changes
      // there's almost certainly a neater way to get these different polarities in here...DJ
      //    ....but this is, I believe, correct!
      int xsign = m == 6 ? -1 : +1;  // endcaps are mirrored in x, not in y

      switch (s) {
        case 0:
          if (collateAlongI)
            stripPos[0] += -xsign * relativePos;
          else
            stripPos[1] += -relativePos;
          break;
        case 1:
          if (collateAlongI)
            stripPos[1] += +relativePos;
          else
            stripPos[0] += -xsign * relativePos;
          break;
        case 2:
          if (collateAlongI)
            stripPos[0] += +xsign * relativePos;
          else
            stripPos[1] += +relativePos;
          break;
        case 3:
          if (collateAlongI)
            stripPos[1] += -relativePos;
          else
            stripPos[0] += +xsign * relativePos;
          break;
      }

    }  // endcap

    //    cout << "new strip pos: " << stripPos[0] << " " << stripPos[1] << " " << stripPos[2] << endl;

    edm4hep::MutableSimCalorimeterHit* newhit = nullptr;

    // check if we have already seen a virtual cell from this strip
    if (newhits.find(cellid) != newhits.end()) {  // already have one from this strip

      edm4hep::MutableSimCalorimeterHit htt = newhits.find(cellid)->second;  // previously made hit in this stirp
      // check that calculated positions are same!
      bool isOK = true;
      for (int i = 0; i < 3; i++) {
        if (fabs(htt.getPosition()[i] - stripPos[i]) > 1) {  // should be identical, but allow a little slop
          isOK = false;
          break;
        }
      }
      if (!isOK) {
        streamlog_out(ERROR) << "strip virtual cell merging: inconsistent position calc!" << std::endl
                             << "please check that the parameter ECAL_strip_nVirtualCells is correctly set..."
                             << getNumberOfVirtualCells() << std::endl
                             << " layer = (k-1) " << km1 << " (gear) " << gearlayer << endl;

        for (int i = 0; i < 3; i++) {
          streamlog_out(DEBUG) << stripPos[i] << " ";
        }
        streamlog_out(DEBUG) << endl;

        for (int i = 0; i < 3; i++) {
          streamlog_out(DEBUG) << htt.getPosition()[i] << " ";
        }
        streamlog_out(DEBUG) << endl;

        // std::cout << "K-1 = " << km1 ;
        // std::cout << "; GEARlay = " << gearlayer;
        // std::cout << "; m = " << m;
        // std::cout << "; s = " << s;
        // std::cout << "; i = " << i_index;
        // std::cout << "; j = " << j_index << std::endl;

        assert(0);
      }

      // set it to the existing one
      newhit = htt;

    } else {  // it's the first hit from this strip
      // create a new hit for the strip
      newhit = &stripcol.create();
      ;
      newhits[cellid] = newhit;
      newHit.setCellID(cellID);
      newhit.setPosition(stripPos);
      newhit.setEnergy(0);  // this is added when we add MC contributions
    }

    // add the MC contriutions
    float eadd(0);
    auto  singleHit = hit.getContributions();

    for (int ij = 0; ij < singleHit.size(); ij++) {
      auto contrib = stripcol.create();
      contrib.setParticle(hit.getParticle(ij));
      contrib.setEnergy(hit.getEnergy(ij) * energyNonuniformityScaling);
      contrib.setTime(hit.getTime(ij));
      contrib.setPDG(hit.getPDG(ij));
      newhit.addToContributions(contrib);
      eadd += hit.getEnergy(ij) * energyNonuniformityScaling;
    }

    float esum(0);
    auto  singlenewHit = newhit.getContributions();
    for (int ij = 0; ij < singlenewHit.size(); ij++) {
      esum += newhit.getEnergy(ij);
    }
  }  // loop over hits

  // move the hits from the temporary storage to the output collection
  for (std::map<int, edm4hep::MutableSimCalorimeterHit*>::iterator itt = newhits.begin(); itt != newhits.end(); itt++) {
    stripcol->addElement(itt->second);
  }

  // return the new collection
  return stripcol;
}

int DDCaloDigi::getNumberOfVirtualCells() {
  // if ( _strip_virt_cells < 0 ) { // not yet set, try to get from gear file
  //   // number of virtual cells in scintillator strip
  //   try {
  //     const gear::GearParameters& pMokka = Global::GEAR->getGearParameters("MokkaParameters");
  //     std::string nVirtualMokkaS = pMokka.getStringVal("Ecal_Sc_number_of_virtual_cells");
  //     _strip_virt_cells = std::atoi( nVirtualMokkaS.c_str() );
  //     streamlog_out (MESSAGE) << "INFO: got number of virtual cells from gear file: " << _strip_virt_cells << endl;
  //   } catch(gear::UnknownParameterException &e) { // not found in gear file, get from steering file parameter...
  //     streamlog_out (MESSAGE) << "INFO: taking number of virtual cells from steering file (not found in gear file): " << _ecalStrip_default_nVirt << endl;
  //     _strip_virt_cells = _ecalStrip_default_nVirt;
  //   }
  // }
  return _strip_virt_cells;
}

std::vector<std::pair<int, int>>& DDCaloDigi::getLayerConfig() {
  // get the layer layout (silicon, scintillator)
  // first element of layerTypes is the preshower

  if (_layerTypes.size() == 0) {
    for (std::string::size_type i = 0; i < _ecalLayout.size(); ++i) {
      // convert each element of string to integer
      // int type = std::atoi( &ccdd ); // this is not well done (must be null-terminated)
      int etype = _ecalLayout[i] - '0';  // this is less obvious, but works...

      switch (etype) {  // these originally defined in Mokka driver SEcalSD04
        case 0:
          _layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          _layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          break;
        case 1:
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          break;
        case 2:
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          break;
        case 3:
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          break;
        case 4:
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          break;
        case 5:
          _layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          break;
        case 6:
          _layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          break;
        case 7:
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          _layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          break;
        case 8:
          _layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          _layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          break;
        default:
          streamlog_out(ERROR) << "ERROR, unknown layer type " << etype << endl;
      }
    }
  }

  return _layerTypes;
}

void DDCaloDigi::checkConsistency(std::string colName, int layer) {
  if (_applyEcalDigi == 0 || _countWarnings > 20)
    return;

  std::pair<int, int> thislayersetup = getLayerProperties(colName, layer);

  if (_applyEcalDigi == 1 && thislayersetup.first != SIECAL) {
    streamlog_out(ERROR) << "collection: " << colName << endl;
    streamlog_out(ERROR) << "you seem to be trying to apply ECAL silicon digitisation to scintillator? Refusing!"
                         << endl;
    streamlog_out(ERROR) << "check setting of ECAL_apply_realistic_digi: " << _applyEcalDigi << endl;
    assert(0);
    _countWarnings++;
  }

  if (_applyEcalDigi == 2 && thislayersetup.first != SCECAL) {
    streamlog_out(ERROR) << "collection: " << colName << endl;
    streamlog_out(ERROR) << "you seem to be trying to apply ECAL scintillator digitisation to silicon? Refusing!"
                         << endl;
    streamlog_out(ERROR) << "check setting of ECAL_apply_realistic_digi: " << _applyEcalDigi << endl;
    assert(0);
    _countWarnings++;
  }

  if (thislayersetup.second != getStripOrientationFromColName(colName)) {
    streamlog_out(ERROR) << "collection: " << colName << endl;
    streamlog_out(ERROR) << "some inconsistency in strip orientation?" << endl;
    streamlog_out(ERROR) << " from collection name: " << getStripOrientationFromColName(colName) << endl;
    streamlog_out(ERROR) << " from layer config string: " << thislayersetup.second << endl;
    _countWarnings++;
  }

  return;
}

std::pair<int, int> DDCaloDigi::getLayerProperties(std::string colName, int layer) {
  std::pair<int, int> thislayersetup(-99, -99);
  std::string         colNameLow(colName);
  std::transform(colNameLow.begin(), colNameLow.end(), colNameLow.begin(), ::tolower);
  if (colNameLow.find("presh") != string::npos) {  // preshower
    if (layer != 0) {
      streamlog_out(WARNING) << "preshower layer with layer index = " << layer << " ??? " << endl;
    } else {
      thislayersetup = getLayerConfig()[layer];
    }
  } else if (colNameLow.find("ring") !=
             string::npos) {  // ecal ring (has no preshower), and is actually always all silicon
    if (layer < int(getLayerConfig().size())) {
      // thislayersetup = getLayerConfig()[layer];
      thislayersetup = std::pair<int, int>(SIECAL, SQUARE);
    } else {
      streamlog_out(WARNING) << "unphysical layer number? " << layer << " " << getLayerConfig().size() << endl;
    }
  } else {  // endcap, barrel
    if (layer + 1 < int(getLayerConfig().size())) {
      thislayersetup = getLayerConfig()[layer + 1];
    } else {
      streamlog_out(WARNING) << "unphysical layer number? " << layer << " " << getLayerConfig().size() << endl;
    }
  }
  return thislayersetup;
}

int DDCaloDigi::getStripOrientationFromColName(std::string colName) {
  int         orientation(-99);
  std::string colNameLow(colName);
  std::transform(colNameLow.begin(), colNameLow.end(), colNameLow.begin(), ::tolower);
  if (colNameLow.find("trans") != string::npos) {
    orientation = STRIP_ALIGN_ACROSS_SLAB;
  } else if (colNameLow.find("long") != string::npos) {
    orientation = STRIP_ALIGN_ALONG_SLAB;
  } else {  // assume square...
    orientation = SQUARE;
    // cout << "WARNING, cannot guess strip orientation! for collection " << colName << endl;
  }
  return orientation;
}
