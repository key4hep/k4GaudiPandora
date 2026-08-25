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
#ifndef PandoraPFAIdeaAlgorithm_h
#define PandoraPFAIdeaAlgorithm_h 1

#include "DualReadoutCaloHitCreator.h"
#include "GeometryCreatorIdea.h"
#include "PfoCreatorIdea.h"
#include "TrackCreatorIdea.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/TrackCollection.h"

#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"

#include "DDRec/DetectorData.h"
#include "Gaudi/Property.h"

#include <string>
#include <vector>

namespace {
class Pandora;
}

// forward declarations for the external clustering algorithm
class ExternalEventParameter;
class ExternalClusterHolder;

struct PandoraPFAIdeaAlgorithm final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection>(
          const edm4hep::TrackCollection&, const std::vector<const edm4hep::CalorimeterHitCollection*>&,
          const std::vector<const edm4hep::ClusterCollection*>&)> {
public:
  PandoraPFAIdeaAlgorithm(const std::string& name, ISvcLocator* svcLoc);
  ~PandoraPFAIdeaAlgorithm() = default;

  StatusCode initialize() override;
  StatusCode finalize() override { return StatusCode::SUCCESS; }

  std::tuple<edm4hep::ClusterCollection, edm4hep::ReconstructedParticleCollection>
  operator()(const edm4hep::TrackCollection& trackColl,
             const std::vector<const edm4hep::CalorimeterHitCollection*>& caloHitColls,
             const std::vector<const edm4hep::ClusterCollection*>& clusterColls) const override;

  const pandora::Pandora* GetPandora() const;

private:
  StatusCode finaliseSteeringParameters();

  SmartIF<IGeoSvc> m_geoSvc;
  ExternalEventParameter* m_extEvtParam = nullptr; ///< external event parameter (pandora::ExternalParameters)
                                                   ///< created by this algo but deleted by Pandora
  std::unique_ptr<ExternalClusterHolder> m_extClusterHolder;

  pandora::Pandora m_pandora;
  std::unique_ptr<GeometryCreatorIdea> m_geometryCreator;
  std::unique_ptr<DualReadoutCaloHitCreator> m_caloHitCreator;
  std::unique_ptr<TrackCreatorIdea> m_trackCreator;
  std::unique_ptr<PfoCreatorIdea> m_pfoCreator;

  GeometryCreatorIdea::Settings m_geometryCreatorSettings;
  DualReadoutCaloHitCreator::Settings m_caloHitCreatorSettings;
  TrackCreatorIdea::Settings m_trackCreatorSettings;
  PfoCreatorIdea::Settings m_pfoCreatorSettings;

  Gaudi::Property<std::string> m_pandoraSettingsXmlFile{this, "PandoraSettingsXmlFile", "",
                                                        "The pandora settings xml file"};

  Gaudi::Property<bool> m_hasHcalEndcap{this, "HasHcalEndcap", true,
      "Whether the geometry has a dual-readout HCAL endcap.  Set false only for a barrel-only "
      "geometry; the calo steering vectors must then drop the endcap subdetector too"};

  // calo hit creator settings
  Gaudi::Property<std::string> m_cherenkovFieldName{this, "CherenkovFieldName", "cherenkov",
                                                    "Name of the cherenkov field in the cellID encoding"};
  // for each detector
  Gaudi::Property<std::vector<uint64_t>> m_systemIDs{
      this, "CaloSystemIDs", {}, "User-given system IDs for the different calorimeters"};
  Gaudi::Property<std::vector<std::string>> m_collectionTypes{
      this,
      "CaloCollectionTypes",
      {},
      "Types of the calo hit collections, either ECAL or HCAL (use ECAL for the monolithic DRC)"};
  Gaudi::Property<std::vector<std::string>> m_layerFieldNames{
      this,
      "CaloLayerFieldNames",
      {},
      "Names of the layer field in the cellID encoding (leave empty if longitudinally unsegmented)"};
  Gaudi::Property<std::vector<std::string>> m_encodingStrings{
      this, "CaloEncodingStrings", {}, "User-given cellID encoding strings (temporary solution)"};
  Gaudi::Property<std::vector<float>> m_cellSizes{
      this, "CaloCellSizes", {}, "Cell sizes in mm (only used for PandoraMonitoring)"};
};

#endif
