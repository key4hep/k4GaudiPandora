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
#ifndef k4GaudiPandora_PfoCreatorIdea_h
#define k4GaudiPandora_PfoCreatorIdea_h 1

#include "PfoCreatorBase.h"

#include "Api/PandoraApi.h"

#include "Gaudi/Algorithm.h"

#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"

/**
 *  @brief  Pfo creator for the IDEA detector: turns the pandora pfos into edm4hep reconstructed
 *          particles and their clusters.
 */
class PfoCreatorIdea : public PfoCreatorBase {
public:
  class Settings {
  public:
    Settings() = default;
    ~Settings() = default;

    /// Dual-readout correction factors, needed to propagate the per-hit energy errors through
    /// E_DR = (S - chi C) / (1 - chi).  Overwritten from the DualReadoutCorrection plugin instance
    /// once pandora has read its settings, so they cannot drift from the chi actually applied; the
    /// defaults only mirror the plugin's own so that a stray unconfigured instance still yields a
    /// sensible scale rather than 1.
    float m_chiEcal = 0.41f;
    float m_chiHcal = 0.31f;
  };

public:
  /**
   *  @brief  Constructor
   *
   *  @param  settings the creator settings
   *  @param  pandora reference to the relevant pandora instance
   *  @param  algorithm reference to the Gaudi::Algorithm to use for message streaming
   */
  PfoCreatorIdea(const Settings& settings, pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm);

  ~PfoCreatorIdea() = default;

  /**
   *  @brief  Create particle flow objects
   *
   *  @param  clusterColl the cluster output collection
   *  @param  aPfoColl the pfo output collection
   */
  pandora::StatusCode CreateParticleFlowObjects(edm4hep::ClusterCollection& clusterColl,
                                                edm4hep::ReconstructedParticleCollection& aPfoColl) const;

private:
  /**
   *  @brief  Create the edm4hep clusters of a pandora pfo and attach them to the reconstructed
   *          particle
   *
   *  @param  pPandoraPfo the address of the pandora pfo
   *  @param  clusterColl the cluster output collection
   *  @param  reconstructedParticle the reconstructed particle to be given the clusters
   */
  pandora::StatusCode AddClustersToRecoParticle(const pandora::ParticleFlowObject* const pPandoraPfo,
                                                edm4hep::ClusterCollection& clusterColl,
                                                edm4hep::MutableReconstructedParticle& reconstructedParticle) const;

  const Settings m_settings; ///< The pfo creator settings
};

#endif // #ifndef k4GaudiPandora_PfoCreatorIdea_h
