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
 *  @file   DDMarlinPandora/include/DDExternalClusteringAlgorithm.h
 *
 *  @brief  Header file for the external clustering algorithm class.
 *
 *  $Log: $
 */
#ifndef DDEXTERNALCLUSTERINGALGORITHM_H
#define DDEXTERNALCLUSTERINGALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <map>

namespace pandora {
class CaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDExternalClusteringAlgorithm class
 */
class DDExternalClusteringAlgorithm : public pandora::Algorithm {
public:
  /**
   *  @brief  Factory class for instantiating algorithm
   */
  class Factory : public pandora::AlgorithmFactory {
  public:
    pandora::Algorithm* CreateAlgorithm() const;
  };

  /**
   *  @brief  Default constructor
   */
  DDExternalClusteringAlgorithm();

private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  typedef std::map<const void*, const pandora::CaloHit*> ParentAddressToCaloHitMap;

  std::string m_externalClusterCollectionName = ""; ///< The collection name for the external clusters
  bool m_flagClustersAsPhotons = true;              ///< Whether to automatically flag new clusters as fixed photons
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm* DDExternalClusteringAlgorithm::Factory::CreateAlgorithm() const {
  return new DDExternalClusteringAlgorithm();
}

#endif // #ifndef DDEXTERNALCLUSTERINGALGORITHM_H
