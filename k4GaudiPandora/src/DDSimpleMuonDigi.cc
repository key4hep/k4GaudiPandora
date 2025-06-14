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

#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include <GaudiKernel/MsgStream.h>

#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/Constants.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include <k4FWCore/MetadataUtils.h>

#include <cstdlib>
#include <functional>
#include <ranges>
#include <utility>

DDSimpleMuonDigi::DDSimpleMuonDigi(const std::string& aName, ISvcLocator* aSvcLoc)
    : MultiTransformer(aName, aSvcLoc,
                       {
                           KeyValues("MUONCollection", {"ECalBarrelCollection"}),
                           KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {KeyValues("MUONOutputCollections", {"CalorimeterHit"}),
                        KeyValues("RelationOutputCollection", {"RelationMuonHit"})}) {}

StatusCode DDSimpleMuonDigi::initialize() {

  m_geoSvc = serviceLocator()->service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to retrieve GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  // Get the number of Layers in the Endcap and Barrel
  size_t layersEndcap = 0, layersBarrel = 0;
  try {
    const auto mainDetector = m_geoSvc->getDetector();
    dd4hep::DetElement theDetector = mainDetector->detector(m_detectorNameBarrel);
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
    const auto mainDetector = m_geoSvc->getDetector();
    dd4hep::DetElement theDetector = mainDetector->detector(m_detectorNameEndcap);
    const dd4hep::rec::LayeredCalorimeterData* yokeEndcapParameters =
        theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();
    layersEndcap = yokeEndcapParameters->layers.size();
  } catch (std::exception& e) {
    debug() << "  oops - no Yoke Endcap available: " << e.what() << std::endl;
  }

  // If the vectors are empty, we are keeping everything
  if (m_layersToKeepBarrelVec.size() > 0) {
    // layers start at 0
    m_useLayersBarrelVec = std::vector<bool>(layersBarrel, false);
    for (const auto k : m_layersToKeepBarrelVec) {
      m_useLayersBarrelVec[k - 1] = true;
    }
  }
  if (m_layersToKeepEndCapVec.size() > 0) {
    // layers start at 0
    m_useLayersEndcapVec = std::vector<bool>(layersEndcap, false);
    for (const auto k : m_layersToKeepEndCapVec) {
      m_useLayersEndcapVec[k - 1] = true;
    }
  }

  const auto collName = inputLocations("MUONCollection")[0];
  const auto encodingString =
      k4FWCore::getParameter<std::string>(collName + "__" + edm4hep::labels::CellIDEncoding, this);
  if (!encodingString) {
    throw std::runtime_error("Encoding string not found for collection: " + collName);
  }
  m_encodingString = encodingString.value();

  return StatusCode::SUCCESS;
}
std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
DDSimpleMuonDigi::operator()(const edm4hep::SimCalorimeterHitCollection& SimCaloHits,
                             const edm4hep::EventHeaderCollection& headers) const {
  debug() << " process event : " << headers[0].getEventNumber() << " - run  " << headers[0].getRunNumber()
          << endmsg; // headers[0].getRunNumber(),headers[0].getEventNumber()

  auto muoncol = edm4hep::CalorimeterHitCollection();
  auto muonRelcol = edm4hep::CaloHitSimCaloHitLinkCollection();

  const auto colName = inputLocations(0)[0];
  CHT::Layout caloLayout = layoutFromString(colName);

  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(m_encodingString);

  for (const auto& hit : SimCaloHits) {
    const auto cellID = hit.getCellID();
    float energy = hit.getEnergy();
    // Get the layer number
    const auto layer = static_cast<size_t>(bitFieldCoder.get(cellID, m_cellIDLayerString));
    // Check if we want to use this later, else go to the next hit
    if (!useLayer(caloLayout, layer))
      continue;
    // Do the digitalization
    float calibr_coeff = 1.;
    calibr_coeff = m_calibrCoeffMuon;
    float hitEnergy = calibr_coeff * energy;
    if (hitEnergy > m_maxHitEnergyMuon) {
      hitEnergy = m_maxHitEnergyMuon;
    }
    if (hitEnergy > m_thresholdMuon) {
      auto calHit = muoncol.create();
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

// If the vectors are empty, we are keeping everything
bool DDSimpleMuonDigi::useLayer(const CHT::Layout caloLayout, const size_t layer) const {
  switch (caloLayout) {
  case CHT::barrel:
    if (layer > m_useLayersBarrelVec.size() || m_useLayersBarrelVec.size() == 0)
      return true;
    return m_useLayersBarrelVec[layer];
  case CHT::endcap:
    if (layer > m_useLayersEndcapVec.size() || m_useLayersEndcapVec.size() == 0)
      return true;
    return m_useLayersEndcapVec[layer];
    // For all other cases, always keep the hit
  default:
    return true;
  }
} // useLayer

float DDSimpleMuonDigi::computeHitTime(const edm4hep::SimCalorimeterHit& h) const {
  // Sort sim hit MC contribution by time.
  // Accumulate the energy from earliest time till the energy
  // threshold is reached. The hit time is then estimated at this position in the array
  using entry_type = std::pair<float, float>;
  std::vector<entry_type> timeToEnergyMapping{};
  timeToEnergyMapping.reserve(h.getContributions().size());

  std::ranges::transform(h.getContributions(), std::back_inserter(timeToEnergyMapping), [](const auto& singleHit) {
    return std::make_pair(singleHit.getTime(), singleHit.getEnergy());
  });
  std::ranges::sort(timeToEnergyMapping, std::less<>{}, [](const auto& entry) { return entry.first; }); // sort by time

  float energySum = 0.f;
  for (const auto& entry : timeToEnergyMapping) {
    energySum += entry.second * m_calibrCoeffMuon;
    if (energySum > m_timeThresholdMuon) {
      return entry.first;
    }
  }

  // default case. That should not happen ...
  return 0.f;
}

DECLARE_COMPONENT(DDSimpleMuonDigi)
