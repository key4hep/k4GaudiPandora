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

#include "edm4hep/ClusterCollection.h"
#include "edm4hep/Cluster.h"

#include "DDExternalClusteringAlgorithm.h"
#include "DDPandoraPFANewProcessor.h"

#include "Pandora/AlgorithmHeaders.h"
#include "XmlHelper.h"

using namespace pandora;

DDExternalClusteringAlgorithm::DDExternalClusteringAlgorithm() : m_flagClustersAsPhotons(true) {}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DDExternalClusteringAlgorithm::Run() {
  try {
    const CaloHitList* pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    if (pCaloHitList->empty())
      return STATUS_CODE_SUCCESS;

    // Get external photon cluster collection
    const edm4hep::EventHeader* const pEventHeader = DDPandoraPFANewProcessor::GetCurrentEvent(&(this->GetPandora()));
    const edm4hep::ClusterCollection* pExternalClusterCollection =
        pEventHeader->getCollection(m_externalClusterCollectionName);
    const unsigned int nExternalClusters(pExternalClusterCollection->size());

    if (0 == nExternalClusters)
      return STATUS_CODE_SUCCESS;

    // Populate pandora parent address to calo hit map
    ParentAddressToCaloHitMap parentAddressToCaloHitMap;

    for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end();
         hitIter != hitIterEnd; ++hitIter) {
      const pandora::CaloHit* const pCaloHit = *hitIter;
      parentAddressToCaloHitMap.insert(ParentAddressToCaloHitMap::value_type(pCaloHit->GetParentAddress(), pCaloHit));
    }

    // Recreate external clusters within the pandora framework
    for (unsigned int iCluster = 0; iCluster < nExternalClusters; ++iCluster) {
      const edm4hep::Cluster* const pExternalCluster = &pExternalClusterCollection->at(iCluster);

      if (nullptr == pExternalCluster)
        throw std::runtime_error("Collection type mismatch");

      const std::vector<edm4hep::ConstCalorimeterHit> calorimeterHitVec = pExternalCluster->getHits();

      const pandora::Cluster* pPandoraCluster = NULL;

      for (auto iter = calorimeterHitVec.begin(), iterEnd = calorimeterHitVec.end(); iter != iterEnd; ++iter) {
        ParentAddressToCaloHitMap::const_iterator pandoraCaloHitIter = parentAddressToCaloHitMap.find(iter->id());

        if (parentAddressToCaloHitMap.end() == pandoraCaloHitIter) {
          continue;
        }

        const pandora::CaloHit* const pPandoraCaloHit = pandoraCaloHitIter->second;

        if (NULL == pPandoraCluster) {
          PandoraContentApi::Cluster::Parameters parameters;
          parameters.m_caloHitList.push_back(pPandoraCaloHit);
          PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                   PandoraContentApi::Cluster::Create(*this, parameters, pPandoraCluster));

          if (m_flagClustersAsPhotons) {
            PandoraContentApi::Cluster::Metadata metadata;
            metadata.m_particleId = PHOTON;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                     PandoraContentApi::Cluster::AlterMetadata(*this, pPandoraCluster, metadata));
          }
        } else {
          PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                   PandoraContentApi::AddToCluster(*this, pPandoraCluster, pPandoraCaloHit));
        }
      }
    }
  } catch (StatusCodeException& statusCodeException) {
    return statusCodeException.GetStatusCode();
  } catch (std::exception& exception) {
    std::cout << "DDExternalClusteringAlgorithm failure: " << exception.what() << std::endl;
    return STATUS_CODE_FAILURE;
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DDExternalClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF(
      STATUS_CODE_SUCCESS, !=,
      XmlHelper::ReadValue(xmlHandle, "ExternalClusterCollectionName", m_externalClusterCollectionName));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "FlagClustersAsPhotons", m_flagClustersAsPhotons));

  return STATUS_CODE_SUCCESS;
}

