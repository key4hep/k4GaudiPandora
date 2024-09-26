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
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/Factories.h"
#include "DDRec/DetectorData.h"
#include "DDRec/MaterialManager.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Constants.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

using namespace dd4hep;
using namespace DDSegmentation;

// Forward Declaration, gets linked in from DDPandoraPFANewProcessor --- FIXME
// now declare here, to be linked from DDPandoraPFANewProcessor
dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag = 0) {
  dd4hep::rec::LayeredCalorimeterData* theExtension = 0;

  dd4hep::Detector&                      mainDetector = dd4hep::Detector::getInstance();
  const std::vector<dd4hep::DetElement>& theDetectors =
      dd4hep::DetectorSelector(mainDetector).detectors(includeFlag, excludeFlag);

  std::cout << " getExtension :  includeFlag: " << dd4hep::DetType(includeFlag)
            << " excludeFlag: " << dd4hep::DetType(excludeFlag) << "  found : " << theDetectors.size()
            << "  - first det: " << theDetectors.at(0).name() << std::endl;

  if (theDetectors.size() != 1) {
    std::stringstream es;
    es << " getExtension: selection is not unique (or empty)  includeFlag: " << dd4hep::DetType(includeFlag)
       << " excludeFlag: " << dd4hep::DetType(excludeFlag) << " --- found detectors : ";
    for (unsigned i = 0, N = theDetectors.size(); i < N; ++i) {
      es << theDetectors.at(i).name() << ", ";
    }
    throw std::runtime_error(es.str());
  }

  theExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

  return theExtension;
};

DDCaloDigi::DDCaloDigi(const std::string& aName, ISvcLocator* aSvcLoc)
    : MultiTransformer(aName, aSvcLoc,
                       {
                           KeyValues("InputCaloHitCollection", {"ECalBarrelCollection"}),
                           KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {KeyValues("OutputCaloHitCollection", {"ECalorimeterHit1"}),
                        KeyValues("RelCollection", {"RelationCaloHit"})}) {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }
  m_geoSvc = serviceLocator()->service("GeoSvc");  // important to initialize m_geoSvc
}

StatusCode DDCaloDigi::initialize() {
  m_countWarnings = 0;

  try {
    dd4hep::rec::LayeredCalorimeterData* ecalBarrelData =
        getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

    dd4hep::rec::LayeredCalorimeterData* ecalEndcapData =
        getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
                     (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& ecalBarrelLayers = ecalBarrelData->layers;
    const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>& ecalEndcapLayers = ecalEndcapData->layers;

    // Determine geometry of ECAL
    int symmetry    = ecalBarrelData->inner_symmetry;
    m_zOfEcalEndcap = ecalEndcapData->extent[2] / dd4hep::mm;

    // Determine ECAL polygon angles
    // Store radial vectors perpendicular to stave layers in m_ecalBarrelStaveDir
    // ASSUMES Mokka Stave numbering 0 = top, then numbering increases anti-clockwise
    if (symmetry > 1) {
      float nFoldSymmetry = static_cast<float>(symmetry);
      float phi0          = ecalBarrelData->phi0 / dd4hep::rad;
      for (int i = 0; i < symmetry; ++i) {
        float phi              = phi0 + i * dd4hep::twopi / nFoldSymmetry;
        m_barrelStaveDir[i][0] = cos(phi);
        m_barrelStaveDir[i][1] = sin(phi);
      }
    }

    for (unsigned int i = 0; i < ecalBarrelLayers.size(); ++i) {
      m_barrelPixelSizeT[i] = ecalBarrelLayers[i].cellSize0;
      m_barrelPixelSizeZ[i] = ecalBarrelLayers[i].cellSize1;
      debug() << "barrel pixel size " << i << " " << m_barrelPixelSizeT[i] << " " << m_barrelPixelSizeZ[i] << std::endl;
    }

    for (unsigned int i = 0; i < ecalEndcapLayers.size(); ++i) {
      m_endcapPixelSizeX[i] = ecalEndcapLayers[i].cellSize0;
      m_endcapPixelSizeY[i] = ecalEndcapLayers[i].cellSize1;
      debug() << "endcap pixel size " << i << " " << m_endcapPixelSizeX[i] << " " << m_endcapPixelSizeY[i] << std::endl;
    }

    m_ecalLayout = m_ecal_deafult_layer_config;
    warning() << "taking layer layout from steering file (FIXME): " << m_ecalLayout << std::endl;

  } catch (std::exception& e) {
    error() << "Could not get ECAL parameters from DD4hep!" << endmsg;
  }

  // Convert ECAL thresholds to GeV units
  if (m_unitThresholdEcal.value().compare("GeV") == 0) {
    // ECAL threshold unit is GeV, do nothing
  } else if (m_unitThresholdEcal.value().compare("MIP") == 0) {
    // ECAL threshold unit is MIP, convert via MIP2GeV
    m_thresholdEcal.value() *= m_calibEcalMip.value();
  } else if (m_unitThresholdEcal.value().compare("px") == 0) {
    // ECAL threshold unit is pixels, convert via MIP2GeV and lightyield
    m_thresholdEcal.value() *= m_ecal_PPD_pe_per_mip.value() * m_calibEcalMip.value();
  } else {
    error() << "Could not identify ECAL threshold unit. Please use \"GeV\", \"MIP\" or \"px\"! Aborting." << endmsg;
    assert(0);
  }

  // Convert HCAL thresholds to GeV units
  if (m_unitThresholdHcal.value().compare("GeV") == 0) {
    // HCAL threshold unit is GeV, do nothing
  } else if (m_unitThresholdHcal.value().compare("MIP") == 0) {
    // HCAL threshold unit is MIP, convert via MIP2GeV
    for (unsigned int i = 0; i < m_thresholdHcal.value().size(); i++) {
      m_thresholdHcal.value()[i] *= m_calibHcalMip.value();
    }
  } else if (m_unitThresholdHcal.value().compare("px") == 0) {
    // HCAL threshold unit is pixels, convert via MIP2GeV and lightyield
    for (unsigned int i = 0; i < m_thresholdHcal.size(); i++) {
      m_thresholdHcal[i] *= m_hcal_PPD_pe_per_mip.value() * m_calibHcalMip.value();
    }
  } else {
    error() << "Could not identify HCAL threshold unit. Please use \"GeV\", \"MIP\" or \"px\"! Aborting." << endmsg;
    assert(0);
  }

  // Set up the scintillator/MPPC digitizer
  m_scEcalDigi = std::unique_ptr<DDScintillatorPpdDigi>(new DDScintillatorPpdDigi());
  m_scEcalDigi->setPEperMIP(m_ecal_PPD_pe_per_mip);
  m_scEcalDigi->setCalibMIP(m_calibEcalMip);
  m_scEcalDigi->setNPix(m_ecal_PPD_n_pixels);
  m_scEcalDigi->setRandomMisCalibNPix(m_ecal_misCalibNpix);
  m_scEcalDigi->setPixSpread(m_ecal_pixSpread);
  m_scEcalDigi->setElecNoise(m_ecal_elec_noise);
  m_scEcalDigi->setElecRange(m_ecalMaxDynMip);
  std::cout << "ECAL sc digi:" << std::endl;
  m_scEcalDigi->printParameters();

  m_scHcalDigi = std::unique_ptr<DDScintillatorPpdDigi>(new DDScintillatorPpdDigi());
  m_scHcalDigi->setPEperMIP(m_hcal_PPD_pe_per_mip);
  m_scHcalDigi->setCalibMIP(m_calibHcalMip);
  m_scHcalDigi->setNPix(m_hcal_PPD_n_pixels);
  m_scHcalDigi->setRandomMisCalibNPix(m_hcal_misCalibNpix);
  m_scHcalDigi->setPixSpread(m_hcal_pixSpread);
  m_scHcalDigi->setElecNoise(m_hcal_elec_noise);
  m_scHcalDigi->setElecRange(m_hcalMaxDynMip);
  std::cout << "HCAL sc digi:" << std::endl;
  m_scHcalDigi->printParameters();

  // Set up the random engines for ECAL and HCAL dead cells: (could use a steering parameter though)
  if (m_deadCellEcal_keep) {
    m_randomEngineDeadCellEcal = new CLHEP::MTwistEngine(0, 0);
  } else {
    m_randomEngineDeadCellEcal = 0;
  }
  if (m_deadCellHcal_keep) {
    m_randomEngineDeadCellHcal = new CLHEP::MTwistEngine(0, 0);
  } else {
    m_randomEngineDeadCellHcal = 0;
  }

  return StatusCode::SUCCESS;
}

retType DDCaloDigi::operator()(const edm4hep::SimCalorimeterHitCollection& simCaloHitCollection,
                               const edm4hep::EventHeaderCollection&       headers) const {
  debug() << " process event : " << headers[0].getEventNumber() << " - run  " << headers[0].getRunNumber() << endmsg;

  auto col = edm4hep::CalorimeterHitCollection();  // create output CalorimeterHit collection
  auto Relcol =
      edm4hep::CaloHitSimCaloHitLinkCollection();  // create relation collection CalorimeterHit-SimCalorimeterHit

  // decide on this event's correlated miscalibration
  //float m_event_correl_miscalib_ecal = CLHEP::RandGauss::shoot(1.0, m_misCalibEcal_correl.value());
  //float m_event_correl_miscalib_hcal = CLHEP::RandGauss::shoot(1.0, m_misCalibHcal_correl.value());

  std::string colName = m_inputLocations[0][0].key();  // take input collection name
  std::cout << "looking for collection: " << colName << std::endl;

  if (colName.find("dummy") != std::string::npos) {
    debug() << "Ignoring input collection name (looks like dummy name)" << colName << endmsg;
  }

  CHT::Layout caloLayout = layoutFromString(colName);

  //FIXME: take this from the metadata !!!
  //std::string initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  std::string initString = "system:5,side:2,module:8,stave:4,layer:9,submodule:4,x:32:-16,y:-16";

  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!
                                                                    // check if decoder contains "layer"

  std::vector<edm4hep::MutableCalorimeterHit*> m_calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS];
  std::vector<int>                             m_calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS];

  int orientation = getStripOrientationFromColName(colName);

  //
  // * Reading Collections of ECAL Simulated Hits *
  //
  if (m_inputColIsECAL) {
    //loop over all SimCalorimeterHits in this collection
    for (const auto& hit : simCaloHitCollection) {
      const int cellID = hit.getCellID();
      float     energy = hit.getEnergy();

      // apply threshold cut
      if (energy > m_thresholdEcal) {
        unsigned int layer  = bitFieldCoder.get(cellID, "layer");
        unsigned int stave  = bitFieldCoder.get(cellID, "stave");
        unsigned int module = bitFieldCoder.get(cellID, "module");

        // check that layer and assumed layer type are compatible
        checkConsistency(colName, layer);

        float x    = hit.getPosition()[0];
        float y    = hit.getPosition()[1];
        float z    = hit.getPosition()[2];
        float r    = sqrt(x * x + y * y + z * z);
        float rxy  = sqrt(x * x + y * y);
        float cost = fabs(z) / r;

        // fill ECAL Layer histograms
        if (z > 0) {
          if (layer == 1)
            ++fEcalLayer1[{x, y}];
          if (layer == 11)
            ++fEcalLayer11[{x, y}];
          if (layer == 21)
            ++fEcalLayer21[{x, y}];
          if (layer == 1)
            ++fEcalRLayer1[rxy];
          if (layer == 11)
            ++fEcalRLayer11[rxy];
          if (layer == 21)
            ++fEcalRLayer21[rxy];
        }

        // evaluate the ECAL calibration coefficient
        float calibr_coeff = 1.;
        if (m_digitalEcal) {
          calibr_coeff = this->digitalEcalCalibCoeff(layer);
          if (m_mapsEcalCorrection) {
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
        // if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff *= m_ecalEndcapCorrectionFactor;
        if (caloLayout != CHT::endcap) {
          calibr_coeff *= m_ecalEndcapCorrectionFactor;  // more robust
        }

        // apply timing cut for ECAL
        if (m_useEcalTiming) {
          float ecalTimeWindowMax;
          if (caloLayout == CHT::barrel) {  // current SimHit is in barrel, use barrel timing cut
            ecalTimeWindowMax = m_ecalBarrelTimeWindowMax;
          } else {  // current simhit is not in barrel, use endcap timing cut
            ecalTimeWindowMax = m_ecalEndcapTimeWindowMax;
          }

          float              dt            = r / 300. - 0.1;
          auto               ecalSingleHit = hit.getContributions();
          int                count         = 0;
          float              eCellInTime   = 0.;
          float              eCellOutput   = 0.;
          const unsigned int n             = ecalSingleHit.size();
          std::vector<bool>  used(n, false);

          for (unsigned int i_t = 0; i_t < n; i_t++) {         // loop over all subhits
            float timei     = ecalSingleHit[i_t].getTime();    // absolute hit timing of current subhit
            float energyi   = ecalSingleHit[i_t].getEnergy();  // energy of current subhit
            float energySum = 0;

            float deltat = 0;
            if (m_ecalCorrectTimesForPropagation) {
              deltat = dt;  //deltat now carries hit timing correction.
            }
            if (timei - deltat > m_ecalTimeWindowMin.value() && timei - deltat < ecalTimeWindowMax) {
              float ecor = energyi * calibr_coeff;
              eCellInTime += ecor;
            }

            if (!used
                    [i_t]) {  //if current subhit has not been merged with previous hits already, take current hit as starting point to merge hits
              // merge with other hits?
              used[i_t] = true;
              for (unsigned int j_t = i_t + 1; j_t < n; j_t++) {  //loop through all hits after current hit
                if (!used[j_t]) {
                  float timej   = ecalSingleHit[j_t].getTime();
                  float energyj = ecalSingleHit[j_t].getEnergy();
                  if (m_ecalSimpleTimingCut) {
                    float deltat_ij = m_ecalCorrectTimesForPropagation ? dt : 0;
                    if (timej - deltat_ij > m_ecalTimeWindowMin && timej - deltat_ij < ecalTimeWindowMax) {
                      energySum += energyj;
                      if (timej < timei) {
                        timei = timej;  //use earliest hit time for simpletimingcut
                      }
                    }
                  } else {
                    float deltat_ij = fabs(timei - timej);
                    //if this subhit is close to current subhit, add this hit's energy to
                    if (deltat_ij < m_ecalDeltaTimeHitResolution) {
                      if (energyj > energyi) {
                        timei = timej;
                      }
                      energyi += energyj;
                      used[j_t] = true;
                    }
                  }
                }
              }
              if (m_ecalSimpleTimingCut) {
                used = std::vector<bool>(n, true);  // mark everything as used to terminate for loop on next run
                energyi += energySum;  // fill energySum back into energyi to have rest of loop behave the same.
              }

              // variables and their behaviour at this point:
              // if SimpleTimingCut == false
              // energyi carries the sum of subhit energies within +- one hcal time resolution - the timecluster energy.
              // timei carries something vaguely similar to the central hit time of the merged subhits

              // if SimpleTimingCut == true
              // energyi carries the sum of subhit energies within timeWindowMin and timeWindowMax
              // timei carries the time of the earliest hit within this window

              if (m_digitalEcal) {
                calibr_coeff = this->digitalEcalCalibCoeff(layer);
                if (m_mapsEcalCorrection) {
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
              if (caloLayout != CHT::barrel) {
                calibr_coeff *= m_ecalEndcapCorrectionFactor;  // more robust
              }

              // fill the ECAL time profile histograms (weighted by the energy)
              fEcal[timei] += energyi * calibr_coeff;
              fEcalC[timei - dt] += energyi * calibr_coeff;
              fEcalC1[timei - dt] += energyi * calibr_coeff;
              fEcalC2[timei - dt] += energyi * calibr_coeff;

              // apply extra energy digitisation effects
              energyi = ecalEnergyDigi(energyi, cellID);  // this only uses the current subhit "timecluster"!

              if (energyi > m_thresholdEcal) {  // now would be the correct time to do threshold comparison
                float timeCor = 0;
                if (m_ecalCorrectTimesForPropagation) {
                  timeCor = dt;
                }
                timei = timei - timeCor;
                if (timei > m_ecalTimeWindowMin &&
                    timei < ecalTimeWindowMax) {  // if current subhit timecluster is within specified timing window,
                                                  // create new CalorimeterHit and add to collections etc.
                  count++;
                  edm4hep::MutableCalorimeterHit calHit = col.create();
                  if (m_ecalGapCorrection != 0) {
                    m_calHitsByStaveLayer[stave][layer].push_back(&calHit);
                    m_calHitsByStaveLayerModule[stave][layer].push_back(module);
                  }
                  calHit.setCellID(cellID);

                  if (m_digitalEcal) {
                    calHit.setEnergy(calibr_coeff);
                  } else {
                    calHit.setEnergy(calibr_coeff * energyi);
                  }

                  eCellOutput += energyi * calibr_coeff;

                  calHit.setTime(timei);
                  calHit.setPosition(hit.getPosition());
                  calHit.setType(CHT(CHT::em, CHT::ecal, caloLayout, layer));

                  // set relation with CaloHitSimCaloHitLinkCollection
                  auto ecaloRel = Relcol.create();
                  ecaloRel.setFrom(calHit);
                  ecaloRel.setTo(hit);

                } else {
                  // if(caloLayout==CHT::barrel) std::cout << " Drop ECAL Barrel hit : " << timei << " " << calibr_coeff*energyi << std::endl;
                }
              }
            }
          }
        } else {  // if timing cut is not used
          edm4hep::MutableCalorimeterHit calHit = col.create();
          if (m_ecalGapCorrection != 0) {
            m_calHitsByStaveLayer[stave][layer].push_back(&calHit);
            m_calHitsByStaveLayerModule[stave][layer].push_back(module);
          }
          float energyi = hit.getEnergy();

          // apply extra energy digitisation effects
          energyi = ecalEnergyDigi(energyi, cellID);

          calHit.setCellID(cellID);
          if (m_digitalEcal) {
            calHit.setEnergy(calibr_coeff);
          } else {
            calHit.setEnergy(calibr_coeff * energyi);
          }
          calHit.setTime(0);
          calHit.setPosition(hit.getPosition());
          calHit.setType(CHT(CHT::em, CHT::ecal, caloLayout, layer));

          // set relation with CaloHitSimCaloHitLinkCollection
          auto ecaloRel = Relcol.create();
          ecaloRel.setFrom(calHit);
          ecaloRel.setTo(hit);
        }  // timing if...else end

        //std::cout << hit->getTimeCont(0) << " count = " << count <<  " E ECAL = " << energyCal << " - " << eCellInTime << " - " << eCellOutput << std::endl;
      }  // energy threshold end
    }    // hits loop end

    // if requested apply gap corrections in ECAL ?
    if (m_ecalGapCorrection != 0) {
      this->fillECALGaps(m_calHitsByStaveLayer, m_calHitsByStaveLayerModule);
    }

    // fill normalisation of HCAL occupancy plots
    for (float x = 15; x < 3000; x += 30) {
      for (float y = 15; y < 3000; y += 30) {
        if (x > 430 || y > 430) {
          float r = sqrt(x * x + y * y);
          fHcalRLayerNorm[r] += 4.;
        }
      }
    }

    // fill normalisation of ECAL occupancy plots
    for (float x = 2.5; x < 3000; x += 5) {
      for (float y = 2.5; y < 3000; y += 5) {
        float r = sqrt(x * x + y * y);
        if (r > 235) {
          fEcalRLayerNorm[r] += 4.;
        }
      }
    }
  }  // end of ECAL digitization

  //
  // * Reading HCAL Collections of Simulated Hits *
  //
  else {
    //loop over all SimCalorimeterHits in this collection
    for (const auto& hit : simCaloHitCollection) {
      float energy = hit.getEnergy();
      //preselect for SimHits with energySum>threshold. Doubtful at least, as lower energy hit might fluctuate up and still be counted
      if (energy > m_thresholdHcal[0] / 2) {
        int          cellID       = hit.getCellID();
        float        calibr_coeff = 1.0;
        unsigned int layer        = bitFieldCoder.get(cellID, "layer");

        // NOTE : for a digital HCAL this does not allow for varying layer thickness
        // with depth - would need a simple mod to make response proportional to layer thickness
        if (m_digitalHcal) {
          calibr_coeff = this->digitalHcalCalibCoeff(caloLayout, energy);
        } else {
          calibr_coeff = this->analogueHcalCalibCoeff(caloLayout, layer);
        }
        // if(fabs(hit->getPosition()[2])>=_zOfEcalEndcap)calibr_coeff*=_hcalEndcapCorrectionFactor;
        if (caloLayout != CHT::endcap)
          calibr_coeff *= m_hcalEndcapCorrectionFactor;  // more robust, is applied to ALL hits outside of barrel.

        //float energyCal = energy*calibr_coeff
        float x = hit.getPosition()[0];
        float y = hit.getPosition()[1];
        float z = hit.getPosition()[2];
        //float r = sqrt(x*x+y*y);
        if (m_useHcalTiming) {
          float hcalTimeWindowMax;
          if (caloLayout == CHT::barrel) {  //current SimHit is in barrel, use barrel timing cut
            hcalTimeWindowMax = m_hcalBarrelTimeWindowMax;
          } else {  //current SimHit is not in barrel, use endcap timing cut
            hcalTimeWindowMax = m_hcalEndcapTimeWindowMax;
          }

          float r = sqrt(
              x * x + y * y +
              z * z);  //this is a crude approximation. assumes initial particle originated at the very center of the detector.
          float              dt            = r / 300 - 0.1;  //magic numbers! ~
          auto               hcalSingleHit = hit.getContributions();
          const unsigned int n             = hcalSingleHit.size();

          std::vector<bool> used(n, false);

          int count = 0;

          for (unsigned int i_t = 0; i_t < n; i_t++) {         // loop over all subhits
            float timei     = hcalSingleHit[i_t].getTime();    //absolute hit timing of current subhit
            float energyi   = hcalSingleHit[i_t].getEnergy();  //energy of current subhit
            float energySum = 0;
            //float deltat = 0; ???
            //if(_hcalCorrectTimesForPropagation)deltat=dt;???
            //deltat now carries hit timing correction. ????
            //std::cout <<"outer:" << i << " " << n << std::endl; ???

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
                  float timej   = hcalSingleHit[j_t].getTime();
                  float energyj = hcalSingleHit[j_t].getEnergy();
                  //              std::cout << " HCAL  deltat_ij : " << deltat_ij << std::endl;
                  if (m_hcalSimpleTimingCut) {
                    float deltat_ij = m_hcalCorrectTimesForPropagation ? dt : 0;
                    if (timej - deltat_ij > m_hcalTimeWindowMin.value() && timej - deltat_ij < hcalTimeWindowMax) {
                      energySum += energyj;
                      if (timej < timei) {
                        timei = timej;  //use earliest hit time for simpletimingcut
                      }
                    }
                  } else {
                    float deltat_ij = fabs(timei - timej);
                    //if this subhit is close to current subhit, add this hit's energy to timecluster
                    if (deltat_ij < m_hcalDeltaTimeHitResolution) {
                      if (energyj > energyi) {
                        //this is probably not what was intended. i guess this should find the largest hit of one timecluster and use its hittime for the cluster, but instead it compares the current hit energy to the sum of already found hit energies
                        timei = timej;
                      }
                      //std::cout << timei << " - " << timej << std::endl;
                      //std::cout << energyi << " - " << energyj << std::endl;
                      energyi += energyj;
                      used[j_t] = true;
                      //std::cout << timei << " " << energyi << std::endl;
                    }
                  }
                }
              }
              if (m_hcalSimpleTimingCut) {
                used = std::vector<bool>(
                    n, true);          //mark everything as used to terminate loop. this is worse than goto. i'm sorry.
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
              energyi = hcalEnergyDigi(energyi, cellID);  //this only uses the current subhit "timecluster"!

              if (energyi > m_thresholdHcal[0]) {  //now would be the correct time to do threshold comparison
                float timeCor = 0;
                if (m_hcalCorrectTimesForPropagation)
                  timeCor = dt;
                timei = timei - timeCor;
                if (timei > m_hcalTimeWindowMin &&
                    timei <
                        hcalTimeWindowMax) {  //if current subhit timecluster is within specified timing window, create new CalorimeterHit and add to collections etc.
                  count++;
                  edm4hep::MutableCalorimeterHit calHit = col.create();
                  calHit.setCellID(cellID);
                  if (m_digitalHcal) {
                    calHit.setEnergy(calibr_coeff);
                    //eCellOutput+= calibr_coeff;
                  } else {
                    calHit.setEnergy(calibr_coeff * energyi);
                    //eCellOutput+= energyi*calibr_coeff;
                  }
                  calHit.setTime(timei);
                  calHit.setPosition(hit.getPosition());
                  calHit.setType(CHT(CHT::had, CHT::hcal, caloLayout, layer));
                  // Set relation with CaloHitSimCaloHitLinkCollection
                  auto hcaloRel = Relcol.create();
                  hcaloRel.setFrom(calHit);
                  hcaloRel.setTo(hit);
                } else {
                  //std::cout << "Drop HCAL hit : " << timei << " " << calibr_coeff*energyi << std::endl; ???
                }
              }
            }
          }       // end loop over all subhits
        } else {  // if timing cuo is not used
          edm4hep::MutableCalorimeterHit calHit = col.create();
          calHit.setCellID(cellID);
          float energyi = hit.getEnergy();

          // apply realistic digitisation
          energyi = hcalEnergyDigi(energyi, cellID);

          if (m_digitalHcal) {
            calHit.setEnergy(calibr_coeff);
          } else {
            calHit.setEnergy(calibr_coeff * energyi);
          }
          //eCellOutput+= energyi*calibr_coeff;
          calHit.setTime(0);
          calHit.setPosition(hit.getPosition());
          calHit.setType(CHT(CHT::had, CHT::hcal, caloLayout, layer));

          auto hcaloRel = Relcol.create();
          hcaloRel.setFrom(calHit);
          hcaloRel.setTo(hit);
        }
        // std::cout << hit->getTimeCont(0) << " count = " << count <<  " EHCAL = " << energyCal << " - " << eCellInTime << " - " << eCellOutput << std::endl;
      }
    }
  }  // end of HCAL digitization

  return std::make_tuple(std::move(col), std::move(Relcol));
}

StatusCode DDCaloDigi::finalize() {
  //delete randomengines if needed
  if (m_randomEngineDeadCellHcal != 0) {
    delete m_randomEngineDeadCellHcal;
  }

  if (m_randomEngineDeadCellEcal != 0) {
    delete m_randomEngineDeadCellEcal;
  }

  return StatusCode::SUCCESS;
}

void DDCaloDigi::fillECALGaps(
    std::vector<edm4hep::MutableCalorimeterHit*> m_calHitsByStaveLayer[MAX_STAVES][MAX_LAYERS],
    std::vector<int>                             m_calHitsByStaveLayerModule[MAX_STAVES][MAX_LAYERS]) const {
  const float slop = 0.25;  // (mm)
  // Loop over hits in the Barrel
  // For each layer calculated differences in hit positions
  // Look for gaps based on expected separation of adjacent hits
  // loop over staves and layers

  for (int is = 0; is < MAX_STAVES; ++is) {
    for (int il = 0; il < MAX_LAYERS; ++il) {
      if (m_calHitsByStaveLayer[is][il].size() > 1) {
        // compare all pairs of hits just once (j>i)

        for (unsigned int i = 0; i < m_calHitsByStaveLayer[is][il].size() - 1; ++i) {
          edm4hep::MutableCalorimeterHit* hiti    = m_calHitsByStaveLayer[is][il][i];
          int                             modulei = m_calHitsByStaveLayerModule[is][il][i];
          float                           xi      = hiti->getPosition()[0];
          float                           yi      = hiti->getPosition()[1];
          float                           zi      = hiti->getPosition()[2];

          for (unsigned int j = i + 1; j < m_calHitsByStaveLayer[is][il].size(); ++j) {
            edm4hep::MutableCalorimeterHit* hitj    = m_calHitsByStaveLayer[is][il][j];
            int                             modulej = m_calHitsByStaveLayerModule[is][il][j];
            float                           xj      = hitj->getPosition()[0];
            float                           yj      = hitj->getPosition()[1];
            float                           zj      = hitj->getPosition()[2];
            float                           dz      = fabs(zi - zj);
            // *** BARREL CORRECTION ***
            if (fabs(zi) < m_zOfEcalEndcap && fabs(zj) < m_zOfEcalEndcap) {
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
              float zminm = 1.0 * m_barrelPixelSizeZ[il] - slop;
              float zmin  = 1.0 * m_barrelPixelSizeZ[il] + slop;
              float zmax  = 2.0 * m_barrelPixelSizeZ[il] - slop;
              float tminm = 1.0 * m_barrelPixelSizeT[il] - slop;
              float tmin  = 1.0 * m_barrelPixelSizeT[il] + slop;
              float tmax  = 2.0 * m_barrelPixelSizeT[il] - slop;

              // criteria for gaps
              // WOULD BE BETTER TO USE GEAR TO CHECK GAPS ARE OF EXPECTED SIZE
              if (dz > zmin && dz < zmax && dt < tminm)
                zgap = true;
              if (dz < zminm && dt > tmin && dt < tmax)
                tgap = true;
              if (dz > zmin && dz < zmax && dt > tmin && dt < tmax)
                ztgap = true;

              if (modulei != modulej) {
                if (dz > zmin && dz < 3.0 * m_barrelPixelSizeZ[il] - slop && dt < tmin)
                  mgap = true;
              }

              // found a gap now apply a correction based on area of gap/area of pixel
              if (zgap || tgap || ztgap || mgap) {
                float ecor = 1.;
                float f    = m_ecalGapCorrectionFactor;  // fudge
                if (mgap)
                  f = m_ecalModuleGapCorrectionFactor;
                if (zgap || mgap)
                  ecor = 1. + f * (dz - m_barrelPixelSizeZ[il]) / 2. / m_barrelPixelSizeZ[il];
                if (tgap)
                  ecor = 1. + f * (dt - m_barrelPixelSizeT[il]) / 2. / m_barrelPixelSizeT[il];
                if (ztgap)
                  ecor = 1. + f * (dt - m_barrelPixelSizeT[il]) * (dz - m_barrelPixelSizeZ[il]) / 4. /
                                  m_barrelPixelSizeT[il] / m_barrelPixelSizeZ[il];
                float ei = hiti->getEnergy() * ecor;
                float ej = hitj->getEnergy() * ecor;
                hiti->setEnergy(ei);
                hitj->setEnergy(ej);
              }

              // *** ENDCAP CORRECTION ***
            } else if (fabs(zi) > m_zOfEcalEndcap && fabs(zj) > m_zOfEcalEndcap && dz < 100) {
              float dx    = fabs(xi - xj);
              float dy    = fabs(yi - yj);
              bool  xgap  = false;
              bool  ygap  = false;
              bool  xygap = false;
              // criteria gaps in the z and t direction

              // x and y need to be swapped in different staves of endcap.
              float pixsizex, pixsizey;
              if (is % 2 == 1) {
                pixsizex = m_endcapPixelSizeY[il];
                pixsizey = m_endcapPixelSizeX[il];
              } else {
                pixsizex = m_endcapPixelSizeX[il];
                pixsizey = m_endcapPixelSizeY[il];
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
                // std::cout <<"NewLDCCaloDigi found endcap gap, adjusting energy! " << xgap << " " << ygap << " " << xygap << " , " << il << std::endl;
                // std::cout << "stave " << is <<  " layer " << il << std::endl;
                // std::cout << "  dx, dy " << dx<< " " << dy << " , sizes = " << pixsizex << " " << pixsizey << std::endl;
                // std::cout << " xmin... " << xmin << " " << xminm << " " << xmax << " ymin... " << ymin << " " << yminm << " " << ymax << std::endl;

                // found a gap make correction
                float ecor = 1.;
                float f    = m_ecalGapCorrectionFactor;  // fudge
                if (xgap)
                  ecor = 1. + f * (dx - pixsizex) / 2. / pixsizex;
                if (ygap)
                  ecor = 1. + f * (dy - pixsizey) / 2. / pixsizey;
                if (xygap)
                  ecor = 1. + f * (dx - pixsizex) * (dy - pixsizey) / 4. / pixsizex / pixsizey;

                // std::cout << "correction factor = " << ecor << std::endl;

                hiti->setEnergy(hiti->getEnergy() * ecor);
                hitj->setEnergy(hitj->getEnergy() * ecor);
              }
            }
          }
        }
      }
    }
  }

  return;
}

float DDCaloDigi::digitalHcalCalibCoeff(CHT::Layout caloLayout, float energy) const {
  float        calib_coeff = 0;
  unsigned int ilevel      = 0;
  for (unsigned int ithresh = 1; ithresh < m_thresholdHcal.size(); ithresh++) {
    // Assume!!!  hit energies are stored as floats, i.e. 1, 2 or 3
    if (energy > m_thresholdHcal[ithresh])
      ilevel = ithresh;  // ilevel = 0 , 1, 2
  }

  switch (caloLayout) {
    case CHT::barrel:
      if (ilevel > m_calibrCoeffHcalBarrel.value().size() - 1) {
        error() << " Semi-digital level " << ilevel << " greater than number of HCAL Calibration Constants ("
                << m_calibrCoeffHcalBarrel.value().size() << ")" << std::endl;
      } else {
        calib_coeff = m_calibrCoeffHcalBarrel.value()[ilevel];
      }
      break;
    case CHT::endcap:
      if (ilevel > m_calibrCoeffHcalEndcap.value().size() - 1) {
        error() << " Semi-digital level " << ilevel << " greater than number of HCAL Calibration Constants ("
                << m_calibrCoeffHcalEndcap.value().size() << ")" << std::endl;
      } else {
        calib_coeff = m_calibrCoeffHcalEndcap.value()[ilevel];
      }
      break;
    case CHT::plug:
      if (ilevel > m_calibrCoeffHcalOther.value().size() - 1) {
        error() << " Semi-digital level " << ilevel << " greater than number of HCAL Calibration Constants ("
                << m_calibrCoeffHcalOther.value().size() << ")" << std::endl;
      } else {
        calib_coeff = m_calibrCoeffHcalOther.value()[ilevel];
      }
      break;
    default:
      error() << " Unknown HCAL Hit Type " << std::endl;
      break;
  }

  return calib_coeff;
}

float DDCaloDigi::analogueHcalCalibCoeff(CHT::Layout caloLayout, int layer) const {
  float calib_coeff = 0;

  for (unsigned int k(0); k < m_hcalLayers.size(); ++k) {
    int min, max;
    if (k == 0)
      min = 0;
    else
      min = m_hcalLayers[k - 1];

    max = m_hcalLayers[k];
    if (layer >= min && layer < max) {
      switch (caloLayout) {
        case CHT::barrel:
          calib_coeff = m_calibrCoeffHcalBarrel[k];
          break;
        case CHT::endcap:
          calib_coeff = m_calibrCoeffHcalEndcap[k];
          break;
        case CHT::plug:
        case CHT::ring:
          calib_coeff = m_calibrCoeffHcalOther[k];
          break;
        default:
          error() << " Unknown HCAL Hit Type " << std::endl;
          break;
      }
    }
  }

  return calib_coeff;
}

float DDCaloDigi::digitalEcalCalibCoeff(int layer) const {
  float calib_coeff = 0.;

  for (unsigned int k(0); k < m_ecalLayers.size(); ++k) {
    int min, max;
    if (k == 0)
      min = 0;
    else
      min = m_ecalLayers[k - 1];

    max = m_ecalLayers[k];
    if (layer >= min && layer < max) {
      calib_coeff = m_calibrCoeffEcal[k];
      break;
    }
  }
  return calib_coeff;
}

float DDCaloDigi::analogueEcalCalibCoeff(int layer) const {
  float calib_coeff = 0;

  // retrieve calibration constants
  for (unsigned int k(0); k < m_ecalLayers.size(); ++k) {
    int min, max;
    if (k == 0) {
      min = 0;
    } else {
      min = m_ecalLayers[k - 1];
    }
    max = m_ecalLayers[k];
    if (layer >= min && layer < max) {
      calib_coeff = m_calibrCoeffEcal[k];
      break;
    }
  }
  return calib_coeff;
}

float DDCaloDigi::ecalEnergyDigi(float energy, int id) const {
  // some extra digi effects (daniel)
  // controlled by m_applyEcalDigi = 0 (none), 1 (silicon), 2 (scintillator)

  // small update for time-constant uncorrelated miscalibrations. DJ, Jan 2015

  float e_out = energy;
  if (m_applyEcalDigi == 1) {
    e_out = siliconDigi(energy);  // silicon digi
  } else if (m_applyEcalDigi == 2) {
    e_out = scintillatorDigi(energy, true);  // scintillator digi
  }
  // add electronics dynamic range
  // Sept 2015: Daniel moved this to the ScintillatorDigi part, so it is applied before unfolding of sipm response
  // if (_ecalMaxDynMip>0) e_out = min (e_out, _ecalMaxDynMip*_calibEcalMip); ???

  // random miscalib
  if (m_misCalibEcal_uncorrel > 0) {
    float miscal(0);
    if (m_misCalibEcal_uncorrel_keep) {
      int id{0};
      if (m_ECAL_cell_miscalibs.find(id) !=
          m_ECAL_cell_miscalibs.end()) {  // this cell was previously seen, and a miscalib stored
        miscal = m_ECAL_cell_miscalibs.at(id);
      } else {  // we haven't seen this one yet, get a miscalib for it
        miscal = CLHEP::RandGauss::shoot(1.0, m_misCalibEcal_uncorrel);
        // FIXME: this is storing miscalibration globally for a run???
        // FIXME _ECAL_cell_miscalibs[id] = miscal; ???
      }
    } else {
      miscal = CLHEP::RandGauss::shoot(1.0, m_misCalibEcal_uncorrel);
    }
    e_out *= miscal;
  }

  if (m_misCalibEcal_correl > 0)
    e_out *= m_event_correl_miscalib_ecal;

  // random cell kill
  if (m_deadCellFractionEcal > 0) {
    if (m_deadCellEcal_keep == true) {
      int id;

      if (m_ECAL_cell_dead.find(id) != m_ECAL_cell_dead.end()) {  // this cell was previously seen
        if (m_ECAL_cell_dead.at(id) == true) {
          e_out = 0;
        }
      } else {  // we haven't seen this one yet, get a miscalib for it
        bool thisDead = (CLHEP::RandFlat::shoot(m_randomEngineDeadCellEcal, .0, 1.0) < m_deadCellFractionEcal);
        // FIXME global map ???
        //_ECAL_cell_dead[id] = thisDead; ???
        if (thisDead == true) {
          e_out = 0;
        }
      }

    } else {
      if (CLHEP::RandFlat::shoot(0.0, 1.0) < m_deadCellFractionEcal)
        e_out = 0;
    }
  }

  return e_out;
}

float DDCaloDigi::hcalEnergyDigi(float energy, int id) const {
  // some extra digi effects (daniel)
  // controlled by _applyHcalDigi = 0 (none), 1 (scintillator/SiPM)

  // small update for time-constant uncorrelated miscalibrations. DJ, Jan 2015

  float e_out(energy);
  if (m_applyHcalDigi == 1)
    e_out = scintillatorDigi(energy, false);  // scintillator digi

  // add electronics dynamic range
  // Sept 2015: Daniel moved this to the ScintillatorDigi part, so it is applied before unfolding of sipm response
  //  if (_hcalMaxDynMip>0) e_out = min (e_out, _hcalMaxDynMip*_calibHcalMip);

  // random miscalib
  //  if (_misCalibHcal_uncorrel>0) e_out*=CLHEP::RandGauss::shoot( 1.0, _misCalibHcal_uncorrel );
  if (m_misCalibHcal_uncorrel > 0) {
    float miscal(0);
    if (m_misCalibHcal_uncorrel_keep) {
      int id;
      if (m_HCAL_cell_miscalibs.find(id) !=
          m_HCAL_cell_miscalibs.end()) {  // this cell was previously seen, and a miscalib stored
        miscal = m_HCAL_cell_miscalibs.at(id);
      } else {  // we haven't seen this one yet, get a miscalib for it
        miscal = CLHEP::RandGauss::shoot(1.0, m_misCalibHcal_uncorrel);
        // FIXME: same as above
        //_HCAL_cell_miscalibs[id] = miscal; ???
      }
    } else {
      miscal = CLHEP::RandGauss::shoot(1.0, m_misCalibHcal_uncorrel);
    }
    e_out *= miscal;
  }

  if (m_misCalibHcal_correl > 0)
    e_out *= m_event_correl_miscalib_hcal;

  // random cell kill
  if (m_deadCellFractionHcal > 0) {
    if (m_deadCellHcal_keep == true) {
      int id;

      if (m_HCAL_cell_dead.find(id) != m_HCAL_cell_dead.end()) {  // this cell was previously seen
        if (m_HCAL_cell_dead.at(id) == true) {
          e_out = 0;
        }
      } else {  // we haven't seen this one yet, get a miscalib for it
        bool thisDead = (CLHEP::RandFlat::shoot(m_randomEngineDeadCellHcal, .0, 1.0) < m_deadCellFractionHcal);
        // FIXME globally dead cell map???
        //FIXME _HCAL_cell_dead[id] = thisDead; ???
        if (thisDead == true) {
          e_out = 0;
        }
      }

    } else {
      if (CLHEP::RandFlat::shoot(0.0, 1.0) < m_deadCellFractionHcal)
        e_out = 0;
    }
  }
  return e_out;
}

float DDCaloDigi::siliconDigi(float energy) const {
  // applies extra digitisation to silicon hits

  // calculate #e-h pairs
  float nehpairs = 1.0e9 * energy / m_ehEnergy;  // check units of energy! m_ehEnergy is in eV, energy in GeV

  // fluctuate it by Poisson
  float smeared_energy = energy * CLHEP::RandPoisson::shoot(nehpairs) / nehpairs;

  // limited electronics dynamic range // Daniel moved electronics dyn range to here
  if (m_ecalMaxDynMip > 0)
    smeared_energy = std::min(smeared_energy, m_ecalMaxDynMip * m_calibEcalMip);

  // add electronics noise
  if (m_ecal_elec_noise > 0)
    smeared_energy += CLHEP::RandGauss::shoot(0, m_ecal_elec_noise * m_calibEcalMip);

  return smeared_energy;
}

float DDCaloDigi::scintillatorDigi(float energy, bool isEcal) const {
  // this applies some extra digitisation to scintillator+PPD hits (PPD=SiPM, MPPC)
  // - poisson fluctuates the number of photo-electrons according to #PEs/MIP
  // - applies PPD saturation according to #pixels
  // Daniel Jeans, Jan/Feb 2014.

  float digiEn(0);
  if (isEcal) {
    digiEn = m_scEcalDigi->getDigitisedEnergy(energy);  //CHECK!!!
  } else {
    digiEn = m_scHcalDigi->getDigitisedEnergy(energy);
  }
  return digiEn;
}

//FIXME should return a reference to somewhere else? or do this every event?
std::vector<std::pair<int, int>> DDCaloDigi::getLayerConfig() const {
  // get the layer layout (silicon, scintillator)
  // first element of layerTypes is the preshower
  std::vector<std::pair<int, int>> m_layerTypes{};
  if (m_layerTypes.size() == 0) {
    for (std::string::size_type i = 0; i < m_ecalLayout.size(); ++i) {
      // convert each element of string to integer
      // int type = std::atoi( &ccdd ); // this is not well done (must be null-terminated)
      int etype = m_ecalLayout[i] - '0';  // this is less obvious, but works...

      switch (etype) {  // these originally defined in Mokka driver SEcalSD04
        case 0:
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          break;
        case 1:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          break;
        case 2:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          break;
        case 3:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          break;
        case 4:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          break;
        case 5:
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          break;
        case 6:
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          break;
        case 7:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ALONG_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          break;
        case 8:
          m_layerTypes.push_back(std::pair<int, int>(SCECAL, STRIP_ALIGN_ACROSS_SLAB));
          m_layerTypes.push_back(std::pair<int, int>(SIECAL, SQUARE));
          break;
        default:
          error() << "ERROR, unknown layer type " << etype << std::endl;
      }
    }
  }

  return m_layerTypes;
}

void DDCaloDigi::checkConsistency(std::string colName, int layer) const {
  if (m_applyEcalDigi == 0 || m_countWarnings > 20)
    return;

  std::pair<int, int> thislayersetup = getLayerProperties(colName, layer);

  if (m_applyEcalDigi == 1 && thislayersetup.first != SIECAL) {
    error() << "collection: " << colName << endmsg;
    error() << "You seem to be trying to apply ECAL silicon digitisation to scintillator? Refusing!" << endmsg;
    error() << "Check setting of ECAL_apply_realistic_digi: " << m_applyEcalDigi << endmsg;
    assert(0);
    m_countWarnings++;
  }

  if (m_applyEcalDigi == 2 && thislayersetup.first != SCECAL) {
    error() << "collection: " << colName << endmsg;
    error() << "You seem to be trying to apply ECAL scintillator digitisation to silicon? Refusing!" << endmsg;
    error() << "Check setting of ECAL_apply_realistic_digi: " << m_applyEcalDigi << endmsg;
    assert(0);
    m_countWarnings++;
  }

  if (thislayersetup.second != getStripOrientationFromColName(colName)) {
    error() << "collection: " << colName << endmsg;
    error() << "Some inconsistency in strip orientation?" << endmsg;
    error() << " from collection name: " << getStripOrientationFromColName(colName) << endmsg;
    error() << " from layer config string: " << thislayersetup.second << endmsg;
    m_countWarnings++;
  }

  return;
}

std::pair<int, int> DDCaloDigi::getLayerProperties(std::string const& colName, int layer) const {
  std::pair<int, int> thislayersetup(-99, -99);
  std::string         colNameLow(colName);
  std::transform(colNameLow.begin(), colNameLow.end(), colNameLow.begin(), ::tolower);
  if (colNameLow.find("presh") != std::string::npos) {  // preshower
    if (layer != 0) {
      warning() << "preshower layer with layer index = " << layer << " ??? " << endmsg;
    } else {
      thislayersetup = getLayerConfig()[layer];
    }
  } else if (colNameLow.find("ring") !=
             std::string::npos) {  // ecal ring (has no preshower), and is actually always all silicon
    if (layer < int(getLayerConfig().size())) {
      // thislayersetup = getLayerConfig()[layer];
      thislayersetup = std::pair<int, int>(SIECAL, SQUARE);
    } else {
      warning() << "unphysical layer number? " << layer << " " << getLayerConfig().size() << endmsg;
    }
  } else {  // endcap, barrel
    if (layer + 1 < int(getLayerConfig().size())) {
      thislayersetup = getLayerConfig()[layer + 1];
    } else {
      warning() << "unphysical layer number? " << layer << " " << getLayerConfig().size() << endmsg;
    }
  }
  return thislayersetup;
}

int DDCaloDigi::getStripOrientationFromColName(std::string const& colName) const {
  int         orientation(-99);
  std::string colNameLow(colName);
  std::transform(colNameLow.begin(), colNameLow.end(), colNameLow.begin(), ::tolower);
  if (colNameLow.find("trans") != std::string::npos) {
    orientation = STRIP_ALIGN_ACROSS_SLAB;
  } else if (colNameLow.find("long") != std::string::npos) {
    orientation = STRIP_ALIGN_ALONG_SLAB;
  } else {  // assume square...
    orientation = SQUARE;
    //std::cout << "WARNING, cannot guess strip orientation! for collection " << colName << std::endl;
  }
  return orientation;
}
