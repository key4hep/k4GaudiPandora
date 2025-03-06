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

#include "DDExternalClusteringAlgorithm.h"
#include "DDPandoraPFANewProcessor.h"

#include "Pandora/AlgorithmHeaders.h"

using namespace pandora;

DDExternalClusteringAlgorithm::DDExternalClusteringAlgorithm() : m_flagClustersAsPhotons(true) {}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DDExternalClusteringAlgorithm::Run() {
 // TODO: This style was not compatible with key4hep and not used for Mucol, so I removed it...
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
