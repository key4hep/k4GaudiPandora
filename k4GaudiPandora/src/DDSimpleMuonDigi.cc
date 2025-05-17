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
#include <cctype>
#include <cstdlib>  // abs
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Constants.h"

DECLARE_COMPONENT(DDSimpleMuonDigi)

DDSimpleMuonDigi::DDSimpleMuonDigi(const std::string& aName, ISvcLocator* aSvcLoc)
    : MultiTransformer(aName, aSvcLoc,
                       {
                           KeyValues("MUONCollection", {"ECalBarrelCollection"}),
                           KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {KeyValues("MUONOutputCollections", {"CalorimeterHit"}),
                        KeyValues("RelationOutputCollection", {"RelationMuonHit"})}) {
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }

  m_geoSvc = serviceLocator()->service("GeoSvc");  // important to initialize m_geoSvc
}

StatusCode DDSimpleMuonDigi::initialize() {
  //Get the number of Layers in the Endcap and Barrel
  int layersEndcap = 0, layersBarrel = 0;
  try {
    const auto                                 mainDetector = m_geoSvc->getDetector();
    dd4hep::DetElement                         theDetector  = mainDetector->detector(m_detectorNameBarrel);
    const dd4hep::rec::LayeredCalorimeterData* yokeBarrelParameters =
        theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();

    layersBarrel = yokeBarrelParameters->layers.size();
    if (!yokeBarrelParameters) {
      error() << "oops - yokeBarrelParameters is a null pointer" << endmsg;
    }

  } catch (std::exception& e) {
    debug() << "  oops - no Yoke Barrel available: " << e.what() << std::endl;
  }
  try {
    const auto                                 mainDetector = m_geoSvc->getDetector();
    dd4hep::DetElement                         theDetector  = mainDetector->detector(m_detectorNameEndcap);
    const dd4hep::rec::LayeredCalorimeterData* yokeEndcapParameters =
        theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();
    layersEndcap = yokeEndcapParameters->layers.size();
  } catch (std::exception& e) {
    debug() << "  oops - no Yoke Endcap available: " << e.what() << std::endl;
  }

  //If the vectors are empty, we are keeping everything
  if (m_layersToKeepBarrelVec.size() > 0) {
    //layers start at 0
    m_useLayersBarrelVec = std::vector<bool>(layersBarrel, false);
    for (int k : m_layersToKeepBarrelVec) {
      m_useLayersBarrelVec[k - 1] = true;
    }
    // for the check
    //for (auto elem : m_useLayersBarrelVec) { std::cout << "m_useLayersBarrelVec " << elem << std::endl; }
  }
  if (m_layersToKeepEndCapVec.size() > 0) {
    //layers start at 0
    m_useLayersEndcapVec = std::vector<bool>(layersEndcap, false);
    for (int k : m_layersToKeepEndCapVec) {
      m_useLayersEndcapVec[k - 1] = true;
    }
    // just for check
    //for (auto elem : m_useLayersEndcapVec) { std::cout << "m_useLayersEndCapVec " << elem << std::endl; }
  }
  return StatusCode::SUCCESS;
}
std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection> DDSimpleMuonDigi::operator()(
    const edm4hep::SimCalorimeterHitCollection& SimCaloHits, const edm4hep::EventHeaderCollection& headers) const {
  debug() << " process event : " << headers[0].getEventNumber() << " - run  " << headers[0].getRunNumber()
          << endmsg;  // headers[0].getRunNumber(),headers[0].getEventNumber()

  auto muoncol    = edm4hep::CalorimeterHitCollection();
  auto muonRelcol = edm4hep::CaloHitSimCaloHitLinkCollection();

  std::string initString;

  std::string colName    = m_muonCollections;
  CHT::Layout caloLayout = layoutFromString(colName);

  //auto col   = headers[0].getCollection(m_muonCollections[i].c_str());
  initString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(initString);  // check!
  // check if decoder contains "layer"

  for (const auto& hit : SimCaloHits) {
    const auto cellID = hit.getCellID();
    float     energy = hit.getEnergy();
    //Get the layer number
    unsigned int layer = bitFieldCoder.get(cellID, "layer");
    //Check if we want to use this later, else go to the next hit
    if (!useLayer(caloLayout, layer))
      continue;
    //Do the digitalization
    float calibr_coeff = 1.;
    calibr_coeff       = m_calibrCoeffMuon;
    float hitEnergy    = calibr_coeff * energy;
    if (hitEnergy > m_maxHitEnergyMuon) {
      hitEnergy = m_maxHitEnergyMuon;
    }
    if (hitEnergy > m_thresholdMuon) {
      edm4hep::MutableCalorimeterHit calHit = muoncol.create();
      calHit.setCellID(cellID);
      calHit.setEnergy(hitEnergy);
      calHit.setPosition(hit.getPosition());
      calHit.setType(CHT(CHT::muon, CHT::yoke, caloLayout, layer));
      calHit.setTime(computeHitTime(hit));
      auto muonRel = muonRelcol.create();
      muonRel.setFrom(calHit);
      muonRel.setTo(hit);
    }
  }

  return std::make_tuple(std::move(muoncol), std::move(muonRelcol));
}

//StatusCode DDSimpleMuonDigi::finalize() { return StatusCode::SUCCESS; }

//If the vectors are empty, we are keeping everything
bool DDSimpleMuonDigi::useLayer(CHT::Layout caloLayout, unsigned int layer) const {
  switch (caloLayout) {
    case CHT::barrel:
      if (layer > m_useLayersBarrelVec.size() || m_useLayersBarrelVec.size() == 0)
        return true;
      return m_useLayersBarrelVec[layer];  //break not needed, because of return
    case CHT::endcap:
      if (layer > m_useLayersEndcapVec.size() || m_useLayersEndcapVec.size() == 0)
        return true;
      return m_useLayersEndcapVec[layer];  //break not needed, because of return
      //For all other cases, always keep the hit
    default:
      return true;
  }
}  //useLayer

float DDSimpleMuonDigi::computeHitTime(const edm4hep::SimCalorimeterHit& h) const {
  // Sort sim hit MC contribution by time.
  // Accumulate the energy from earliest time till the energy
  // threshold is reached. The hit time is then estimated at this position in the array
  using entry_type = std::pair<float, float>;
  std::vector<entry_type> timeToEnergyMapping{};
  //for(const auto& ih : h) {
  auto singleHit = h.getContributions();

  const unsigned int nContribs = singleHit.size();
  timeToEnergyMapping.reserve(nContribs);

  for (unsigned int c = 0; c < nContribs; ++c) {
    timeToEnergyMapping.push_back({singleHit[c].getTime(), singleHit[c].getEnergy()});
  }
  std::sort(timeToEnergyMapping.begin(), timeToEnergyMapping.end(),
            [this](entry_type& lhs, entry_type& rhs) { return lhs.first < rhs.first; });
  float energySum = 0.f;
  for (auto& entry : timeToEnergyMapping) {
    energySum += entry.second * m_calibrCoeffMuon;
    if (energySum > m_timeThresholdMuon) {
      return entry.first;
    }
  }
  //}
  // default case. That should not happen ...
  return 0.f;
}
