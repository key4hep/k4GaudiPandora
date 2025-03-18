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
#ifndef DDSCINTILLATORPPDDIGI_H
#define DDSCINTILLATORPPDDIGI_H

#include "CLHEP/Random/MTwistEngine.h"

#include <memory>

//GaudiKernel::MsgStream
class MsgStream;

class DDScintillatorPpdDigi {
public:
  DDScintillatorPpdDigi();
  ~DDScintillatorPpdDigi() {}

  // expected # photoelectrons / MIP
  void setPEperMIP(float x) { m_PEperMIP = x; }

  // calibration factor from input hit energy to MIPs
  void setCalibMIP(float x) { m_calibMIP = x; }

  // # pixels of PPD
  void setNPix(int x) { m_Npix = x; }

  // random miscalibration of total # pixels (as a fraction of pixel number: 0.05 = 5% miscalibration)
  void setRandomMisCalibNPix(float x) { m_misCalibNpix = x; }

  // spread in pixel capacitance (as a fraction: 0.05 = 5% spread)
  void setPixSpread(float x) { m_pixSpread = x; }

  // electronics noise (in MIP units)
  void setElecNoise(float x) { m_elecNoise = x; }

  // electronics dynamic range (in MIP units)
  void setElecRange(float x) { m_elecMaxDynRange_MIP = x; }

  float getDigitisedEnergy(float energy, std::unique_ptr<CLHEP::MTwistEngine>&);

  void printParameters(MsgStream&);

private:
  float m_PEperMIP            = -99;
  float m_calibMIP            = -99;
  float m_Npix                = -99;
  float m_misCalibNpix        = 0;
  float m_pixSpread           = 0;
  float m_elecNoise           = 0;
  float m_elecMaxDynRange_MIP = 0;
};

#endif
