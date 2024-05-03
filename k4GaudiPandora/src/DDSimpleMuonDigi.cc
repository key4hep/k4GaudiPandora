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
#include "DDSimpleMuonDigi.h"
#include "GaudiKernel/MsgStream.h"
#include "k4FWCore/BaseClass.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Constants.h"
#include <cctype>
#include <cstdlib>  // abs

DECLARE_COMPONENT(DDSimpleMuonDigi)

DDSimpleMuonDigi::DDSimpleMuonDigi(const std::string& aName, ISvcLocator* aSvcLoc)
: MultiTransformer(aName, aSvcLoc,
		   {
		     KeyValues("MUONCollections", {"SimCalorimeterHit"}) KeyValues("HeaderName", {"EventHeader"}),
		   },
		   {KeyValues("MUONOutputCollections", {"CalorimeterHit"})
		    KeyValues("RelationOutputCollection", {"RelationMuonHit"})}){

  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }


}

          
StatusCode DDSimpleMuonDigi::initialize() {
  //Get the number of Layers in the Endcap
  int layersEndcap = 0, layersBarrel = 0;
  try {
    const auto detector = m_geoSvc->getDetector();
    //dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();
    dd4hep::DetElement theDetector = mainDetector.detector(m_detectorNameBarrel);
    const dd4hep::rec::LayeredCalorimeterData* yokeBarrelParameters =
      theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();
    layersBarrel = yokeBarrelParameters->layers.size();
  } catch (std::exception& e) {
    streamlog_out(WARNING) << "  oops - no Yoke Barrel available: " << e.what() << std::endl;
  }
  try {
    const auto detector = m_geoSvc->getDetector();
    //dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();
    dd4hep::DetElement  theDetector = mainDetector.detector(m_detectorNameEndcap);
    const dd4hep::rec::LayeredCalorimeterData* yokeEndcapParameters =
      theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();
    layersEndcap = yokeEndcapParameters->layers.size();
  } catch (std::exception& e) {
    streamlog_out(WARNING) << "  oops - no Yoke Endcap available: " << e.what() << std::endl;
  }

  //If the vectors are empty, we are keeping everything
  if (layersToKeepBarrelVec.size() > 0) {
    //layers start at 0
    for (int i = 0; i < layersBarrel; ++i) {
      useLayersBarrelVec.push_back(false);
      for (IntVec::iterator iter = layersToKeepBarrelVec.begin(); iter < layersToKeepBarrelVec.end(); ++iter) {
        if (i == *iter - 1) {
          useLayersBarrelVec[i] = true;
          break;
        }
      }
    }
  }

  if (layersToKeepEndcapVec.size() > 0) {
    //layers start at 0
    for (int i = 0; i < layersEndcap; ++i) {
      useLayersEndcapVec.push_back(false);
      for (IntVec::iterator iter = layersToKeepEndcapVec.begin(); iter < layersToKeepEndcapVec.end(); ++iter) {
        if (i == *iter - 1) {
          useLayersEndcapVec[i] = true;
          break;
        }
      }
    }
  }
}
DDSimpleMuonDigi::operator()(const edm4hep::SimCalorimeterHit &  SimCaloHits,
                             const edm4hep::EventHeaderCollection & headers) const {
  debug() << " process event : " << headers[0].getEventNumber() << " - run  " << headers[0].getRunNumber()
          << endmsg;  // headers[0].getRunNumber(),headers[0].getEventNumber()

  auto muoncol  = edm4hep::CalorimeterHitCollection();
  auto chschcol = edm4hep::MCRecoCaloAssociationCollection();
  auto muonRelcol = chschcol.create();

  std::string initString;
  for (unsigned int i(0); i < muon_Collections.size(); ++i) {
    std::string colName    = muon_Collections[i];
    CHT::Layout caloLayout = layoutFromString(colName);
    try {
      auto col   = headers[0].getCollection(muon_Collections[i].c_str());
      initString  = m_geoSvc->constantAsString(m_encodingStringVariable.value());
      dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!
      // check if decoder contains "layer"
      std::vector<std::string> fields;
      int                      numElements = col.getNumberOfElements();
    }
    for (const auto& hit : simCaloHits) {
      const int    cellID = hit.getCellID();
      float        energy = hit.getEnergy();
      unsigned int layer  = bitFieldCoder(cellID, _cellIDLayerString);
      if (!useLayer(caloLayout, layer))
        continue;
      float calibr_coeff(1.);
      calibr_coeff    = calibrCoeffMuon;
      float hitEnergy = calibr_coeff * energy;
      if (hitEnergy > maxHitEnergyMuon)
        hitEnergy = maxHitEnergyMuon;
      if (hitEnergy > thresholdMuon) {
        edm4hep::CalorimeterHit& calhit = new CalorimeterHit();
        calhit.setCellID(cellid);
        calhit.setEnergy(hitEnergy);
        calhit.setPosition(hit.getPosition());
        calhit.setType(CHT(CHT::muon, CHT::yoke, caloLayout, idDecoder(hit)[cellIDLayerString]));
        calhit.setTime(computeHitTime(hit));
        calhit.setRawHit(hit);
        muoncol.addElement(calhit);
        muonRelcol.setRec(calhit);
        muonRelcol.setSim(hit);
      }
    }

    catch (std::exception& e) {
    }
  }
  return std::make_tuple(std::move(muoncol), std::move(muonRelcol));
  muoncol.parameters().setValue(edm4hep::CellIDEncoding, initString);
 
}

StatusCode DDSimpleMuonDigi::finalize() {  
  return StatusCode::SUCCESS;

} //fix

bool DDSimpleMuonDigi::useLayer(CHT::Layout caloLayout, unsigned int layer) {
  switch (caloLayout) {
  case CHT::barrel:
    if (layer > useLayersBarrelVec.size() || useLayersBarrelVec.size() == 0)
      return true;
    return useLayersBarrelVec[layer];  //break not needed, because of return
  case CHT::endcap:
    if (layer > useLayersEndcapVec.size() || useLayersEndcapVec.size() == 0)
      return true;
    return useLayersEndcapVec[layer];  //break not needed, because of return
    //For all other cases, always keep the hit
  default:
    return true;
  }
}  //useLayer

float DDSimpleMuonDigi::computeHitTime(const edmp4hep::SimCalorimeterHit* h) const {
  if (nullptr == h) {
    return 0.f;
  }
  // Sort sim hit MC contribution by time.
  // Accumulate the energy from earliest time till the energy
  // threshold is reached. The hit time is then estimated at this position in the array
  using entry_type = std::pair<float, float>;
  std::vector<entry_type> timeToEnergyMapping{};

  const unsigned int nContribs = h->getNMCContributions();
  timeToEnergyMapping.reserve(nContribs);

  for (unsigned int c = 0; c < nContribs; ++c) {
    timeToEnergyMapping.push_back({h->getTimeCont(c), h->getEnergyCont(c)});
  }
  std::sort(timeToEnergyMapping.begin(), timeToEnergyMapping.end(),
            [this](entry_type& lhs, entry_type& rhs) { return lhs.first < rhs.first; });
  float energySum = 0.f;
  for (auto& entry : timeToEnergyMapping) {
    energySum += entry.second * calibrCoeffMuon;
    if (energySum > m_timeThresholdMuon) {
      return entry.first;
    }
  }
  // default case. That should not happen ...
  return 0.f;
}
