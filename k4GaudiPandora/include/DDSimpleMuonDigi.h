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
#ifndef DDSimpleMuonDigi_H
#define DDSimpleMuonDigi_H

#include "CalorimeterHitType.h"

#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/EventHeaderCollection.h>
#include <edm4hep/SimCalorimeterHitCollection.h>

#include <Gaudi/Property.h>

#include <k4FWCore/Transformer.h>
#include <k4Interface/IGeoSvc.h>

#include <string>
#include <tuple>
#include <vector>

struct DDSimpleMuonDigi final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>(
              const edm4hep::SimCalorimeterHitCollection&, const edm4hep::EventHeaderCollection&)> {
  DDSimpleMuonDigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::CaloHitSimCaloHitLinkCollection>
  operator()(const edm4hep::SimCalorimeterHitCollection& simCaloHits,
             const edm4hep::EventHeaderCollection& headers) const override;

private:
  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", "VXD", "Name of the subdetector"};
  Gaudi::Property<std::vector<size_t>> m_layersToKeepBarrelVec{
      this, "KeepBarrelLayersVec", {}, "Vector of Barrel layers to be kept. Layers start at 1!"};
  Gaudi::Property<std::vector<size_t>> m_layersToKeepEndCapVec{
      this, "KeepEndcapLayersVec", {}, "Vector of Endcap layers to be kept. Layers start at 1!"};
  // Gaudi::Property<std::vector<bool>> useLayersBarrelVec{this, "useBarrelLayerVector", false, "whether to use the
  // endcap layer vector"}; Gaudi::Property<std::vector<bool>> useLayersEndcapVec{this, "useEndCapLayerVector", false,
  // "whether to use the EndCap layer vector"};
  Gaudi::Property<std::string> m_cellIDLayerString{this, "CellIDLayerString", "layer",
                                                   "Name of the part of the cellID that holds the layer"};
  Gaudi::Property<float> m_thresholdMuon{this, "MuonThreshold", {0.025f}, "Threshold for muon"};
  Gaudi::Property<float> m_timeThresholdMuon{this, "timethresholdMuon", {0.025f}, "time threshold for muons"};
  Gaudi::Property<float> m_calibrCoeffMuon{this, "CalibrMUON", {120000.0}, "Callibration coefficient of muons"};
  Gaudi::Property<float> m_maxHitEnergyMuon{this, "MaxHitEnergyMUON", {2.0}, "Threshold for maximum muon hit energy"};
  Gaudi::Property<std::string> m_detectorNameBarrel{this, "detectornameB", "YokeBarrel", "Name of the subdetector"};
  Gaudi::Property<std::string> m_detectorNameEndcap{this, "detectornameE", "YokeEndcap",
                                                    "Name of the second subdetector"};

  std::string m_collName{};
  std::string m_encodingString{};
  std::vector<bool> m_useLayersBarrelVec{}, m_useLayersEndcapVec{};
  SmartIF<IGeoSvc> m_geoSvc;

  bool useLayer(const CHT::Layout caloLayout, const size_t layer) const;
  float computeHitTime(const edm4hep::SimCalorimeterHit& h) const;
};
#endif
