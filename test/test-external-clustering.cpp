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

#include "DDExternalClusteringAlgorithm.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"

#include "Api/PandoraApi.h"
#include "Pandora/AlgorithmHeaders.h"
#include "Plugins/PseudoLayerPlugin.h"

#include <iostream>

namespace {
class DummyPseudoLayer : public pandora::PseudoLayerPlugin {
  unsigned int GetPseudoLayer(const pandora::CartesianVector&) const override { return 0; }
  unsigned int GetPseudoLayerAtIp() const override { return 0; }
  pandora::StatusCode ReadSettings(pandora::TiXmlHandle) override { return pandora::STATUS_CODE_SUCCESS; }
};

class CheckCluster : public pandora::Algorithm {
public:
  CheckCluster(const edm4hep::CalorimeterHit& expected, bool& matched) : m_expected(expected), m_matched(matched) {}

  class Factory : public pandora::AlgorithmFactory {
  public:
    Factory(const edm4hep::CalorimeterHit& expected, bool& matched) : m_expected(expected), m_matched(matched) {}
    pandora::Algorithm* CreateAlgorithm() const override { return new CheckCluster(m_expected, m_matched); }

  private:
    const edm4hep::CalorimeterHit& m_expected;
    bool& m_matched;
  };

private:
  pandora::StatusCode ReadSettings(pandora::TiXmlHandle) override { return pandora::STATUS_CODE_SUCCESS; }

  pandora::StatusCode Run() override {
    const pandora::ClusterList* clusters = nullptr;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, clusters));
    if (clusters->size() == 1) {
      pandora::CaloHitList hits;
      clusters->front()->GetOrderedCaloHitList().FillCaloHitList(hits);
      m_matched = hits.size() == 1 && hits.front()->GetParentAddress() == &m_expected;
    }
    return pandora::STATUS_CODE_SUCCESS;
  }

  const edm4hep::CalorimeterHit& m_expected;
  bool& m_matched;
};

pandora::StatusCode createPandoraHit(const pandora::Pandora& pandora, const edm4hep::CalorimeterHit& parent) {
  PandoraApi::CaloHit::Parameters hit;
  hit.m_positionVector = pandora::CartesianVector(100.f, 0.f, 0.f);
  hit.m_expectedDirection = pandora::CartesianVector(1.f, 0.f, 0.f);
  hit.m_cellNormalVector = pandora::CartesianVector(1.f, 0.f, 0.f);
  hit.m_cellGeometry = pandora::RECTANGULAR;
  hit.m_cellSize0 = hit.m_cellSize1 = hit.m_cellThickness = 1.f;
  hit.m_nCellRadiationLengths = hit.m_nCellInteractionLengths = 1.f;
  hit.m_time = 0.f;
  hit.m_inputEnergy = hit.m_mipEquivalentEnergy = hit.m_electromagneticEnergy = hit.m_hadronicEnergy = 1.f;
  hit.m_isDigital = false;
  hit.m_hitType = pandora::ECAL;
  hit.m_hitRegion = pandora::BARREL;
  hit.m_layer = 0u;
  hit.m_isInOuterSamplingLayer = false;
  hit.m_pParentAddress = &parent;
  return PandoraApi::CaloHit::Create(pandora, hit);
}
} // namespace

int main(int argc, char* argv[]) try {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " pandora-settings.xml\n";
    return 1;
  }

  // Before the collections are put in a frame, their first hits share an object ID.
  edm4hep::CalorimeterHitCollection firstCollection, secondCollection;
  const edm4hep::CalorimeterHit first = firstCollection.create();
  const edm4hep::CalorimeterHit second = secondCollection.create();
  if (first.id() != second.id()) {
    std::cerr << "The test requires distinct hits with the same object ID\n";
    return 1;
  }

  edm4hep::ClusterCollection clusters;
  auto cluster = clusters.create();
  cluster.addToHits(second);
  std::vector<std::vector<edm4hep::Cluster>> externalClusters{{cluster}};
  ExternalClusterHolder holder;
  holder.setExternalClusters(&externalClusters);

  bool matched = false;
  pandora::Pandora pandora;
  auto* parameters = new ExternalEventParameter;
  parameters->m_externalClusterHolder = &holder;
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                          PandoraApi::SetExternalParameters(pandora, "DDExternalClustering", parameters));
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                          PandoraApi::RegisterAlgorithmFactory(pandora, "DDExternalClustering",
                                                               new DDExternalClusteringAlgorithm::Factory));
  PANDORA_THROW_RESULT_IF(
      pandora::STATUS_CODE_SUCCESS, !=,
      PandoraApi::RegisterAlgorithmFactory(pandora, "CheckCluster", new CheckCluster::Factory(second, matched)));
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                          PandoraApi::SetPseudoLayerPlugin(pandora, new DummyPseudoLayer));
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(pandora, argv[1]));
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, createPandoraHit(pandora, first));
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, createPandoraHit(pandora, second));
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(pandora));
  if (!matched) {
    std::cerr << "Expected one external cluster containing only the hit from the second collection\n";
    return 1;
  }
  return 0;
} catch (const pandora::StatusCodeException& error) {
  std::cerr << error.ToString() << '\n';
  return 1;
}
