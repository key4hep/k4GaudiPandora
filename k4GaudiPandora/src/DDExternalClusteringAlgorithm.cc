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

/**
 *  @file   DDMarlinPandora/src/DDExternalClusteringAlgorithm.cc
 *
 *  @brief  Implementation of the external clustering algorithm class.
 *
 *  $Log: $
 */

#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"

#include "Pandora/AlgorithmHeaders.h"

#include "DDExternalClusteringAlgorithm.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/AnyDataWrapper.h"
#include "GaudiKernel/IDataProviderSvc.h"

// setters and getters for the external cluster holder
void ExternalClusterHolder::setExternalClusters(std::vector<std::vector<edm4hep::Cluster>>* externalClusters) {
  m_externalClusters = externalClusters;
}

const std::vector<std::vector<edm4hep::Cluster>>& ExternalClusterHolder::getExternalClusters() const {
  return *m_externalClusters;
}

DDExternalClusteringAlgorithm::DDExternalClusteringAlgorithm() : m_flagClustersAsPhotons(false) {}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDExternalClusteringAlgorithm::Run() {
  try {
    // Get Pandora calo hit list
    const pandora::CaloHitList* pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    if (pCaloHitList->empty())
      return pandora::STATUS_CODE_SUCCESS;

    // Populate pandora parent address to calo hit map
    ExternalToPandoraCaloHitMap caloHitMap;

    for (const auto* pCaloHit : *pCaloHitList) {
      const edm4hep::CalorimeterHit* edmCaloHit =
          static_cast<const edm4hep::CalorimeterHit*>(pCaloHit->GetParentAddress());

      caloHitMap.emplace(edmCaloHit->id(), pCaloHit);
    }

    // Recreate external clusters within the pandora framework
    const pandora::ClusterList* pClusterList = nullptr;

    std::string clusterListNameTmp = "ExternalClustersTmp";
    std::string clusterListNameFinal = "ExternalClusters";

    PANDORA_RETURN_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListNameTmp));

    // loop over external cluster collections
    for (const auto& clusterColl : m_externalClusterHolder->getExternalClusters()) {
      // no cluster in this event
      if (clusterColl.empty())
        continue;

      for (const auto& externalCluster : clusterColl) {
        const auto& calorimeterHitVec = externalCluster.getHits();
        pandora::CaloHitList pandoraHitList;

        // find corresponding pandora hits
        // and fill pandora calo hit list
        for (const auto& edmHit : calorimeterHitVec) {
          auto itr = caloHitMap.find(edmHit.id());

          if (itr == caloHitMap.end())
            continue;

          pandoraHitList.push_back(itr->second);
        }

        if (pandoraHitList.empty())
          continue;

        object_creation::ClusterParameters clusterParameters;
        clusterParameters.m_caloHitList = pandoraHitList;

        // create a pandora cluster
        const pandora::Cluster* pPandoraCluster = nullptr;
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                PandoraContentApi::Cluster::Create(*this, clusterParameters, pPandoraCluster));

        if (m_flagClustersAsPhotons) {
          PandoraContentApi::Cluster::Metadata metadata;
          metadata.m_particleId = pandora::PHOTON;
          PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                   PandoraContentApi::Cluster::AlterMetadata(*this, pPandoraCluster, metadata));
        }
      } // loop over external clusters
    } // loop over external cluster collections

    // need these to store the pandora cluster list
    if (!pClusterList->empty()) {
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                               PandoraContentApi::SaveList<pandora::Cluster>(*this, clusterListNameFinal));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                               PandoraContentApi::ReplaceCurrentList<pandora::Cluster>(*this, clusterListNameFinal));
    }
  } catch (pandora::StatusCodeException& statusCodeException) {
    return statusCodeException.GetStatusCode();
  } catch (std::exception& exception) {
    std::cerr << "DDExternalClusteringAlgorithm failure: " << exception.what() << std::endl;
    return pandora::STATUS_CODE_FAILURE;
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDExternalClusteringAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "FlagClustersAsPhotons", m_flagClustersAsPhotons));

  // Get the pointer to the external cluster holder
  // Sanghyun: this is a hacky way that I do this in ReadSettings
  // but GetExternalParameters doesn't allow me to access it more that once
  const ExternalEventParameter* pExternalEventParameter =
      dynamic_cast<const ExternalEventParameter*>(this->GetExternalParameters());

  if (!pExternalEventParameter || !pExternalEventParameter->m_externalClusterHolder) {
    throw std::runtime_error("DDExternalClusteringAlgorithm: External cluster holder not set");
  }

  m_externalClusterHolder = pExternalEventParameter->m_externalClusterHolder;

  return pandora::STATUS_CODE_SUCCESS;
}
