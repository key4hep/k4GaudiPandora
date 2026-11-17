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
#ifndef k4GaudiPandora_PfoCreatorBase_h
#define k4GaudiPandora_PfoCreatorBase_h 1

#include "Api/PandoraApi.h"

#include <Gaudi/Algorithm.h>

namespace edm4hep {
class MutableReconstructedParticle;
}

/**
 *  @brief  Common context and pfo-to-edm4hep transcription shared by all pfo creators.
 *
 *  Implementation sharing only: there are no virtual functions and the destructor is protected, so
 *  a PfoCreatorBase is never held or deleted polymorphically.  Everything that depends on
 *  calorimeter hit semantics lives in the derived classes.
 */
class PfoCreatorBase {
public:
  PfoCreatorBase(const PfoCreatorBase&) = delete;
  PfoCreatorBase& operator=(const PfoCreatorBase&) = delete;
  PfoCreatorBase(PfoCreatorBase&&) = delete;
  PfoCreatorBase& operator=(PfoCreatorBase&&) = delete;

protected:
  /**
   *  @brief  Constructor
   *
   *  @param  pandora reference to the relevant pandora instance
   *  @param  algorithm reference to the Gaudi::Algorithm to use for message streaming
   */
  PfoCreatorBase(pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm);

  ~PfoCreatorBase() = default;

  /**
   *  @brief  Add the tracks of the pandora pfo to the reconstructed particle
   *
   *  @param  pPandoraPfo the address of the pandora pfo
   *  @param  pReconstructedParticle the reconstructed particle to be given the tracks
   */
  void AddTracksToRecoParticle(const pandora::ParticleFlowObject* const pPandoraPfo,
                               edm4hep::MutableReconstructedParticle& pReconstructedParticle) const;

  /**
   *  @brief  Copy momentum, energy, mass, charge and pdg from the pandora pfo
   *
   *  @param  pPandoraPfo the address of the pandora pfo
   *  @param  pReconstructedParticle the reconstructed particle to be given the properties
   */
  void SetRecoParticlePropertiesFromPFO(const pandora::ParticleFlowObject* const pPandoraPfo,
                                        edm4hep::MutableReconstructedParticle& pReconstructedParticle) const;

  pandora::Pandora& m_pandora;         ///< Reference to the pandora object from which to extract the pfos
  const Gaudi::Algorithm& m_algorithm; ///< Reference to the Gaudi algorithm for message streaming
};

#endif // #ifndef k4GaudiPandora_PfoCreatorBase_h
