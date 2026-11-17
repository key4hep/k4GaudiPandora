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

#include "PfoCreatorBase.h"

#include <edm4hep/MutableReconstructedParticle.h>
#include <edm4hep/Track.h>

#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

PfoCreatorBase::PfoCreatorBase(pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm)
    : m_pandora(pandora), m_algorithm(*algorithm) {}

void PfoCreatorBase::AddTracksToRecoParticle(const pandora::ParticleFlowObject* const pPandoraPfo,
                                             edm4hep::MutableReconstructedParticle& pReconstructedParticle) const {
  // Every pandora track was built from an edm4hep track whose address the track creator stored as
  // the parent address, so the association is recovered by casting that address back rather than
  // by searching the input collection.
  for (const auto* pTrack : pPandoraPfo->GetTrackList()) {
    const auto& pLcioTrack = *static_cast<const edm4hep::Track*>(pTrack->GetParentAddress());
    pReconstructedParticle.addToTracks(pLcioTrack);
  }
}

void PfoCreatorBase::SetRecoParticlePropertiesFromPFO(
    const pandora::ParticleFlowObject* const pPandoraPfo,
    edm4hep::MutableReconstructedParticle& reconstructedParticle) const {
  // Straight transcription of the pfo-level quantities.  Nothing is recomputed here: the energy,
  // the momentum and the particle hypothesis were all decided by the pfo construction algorithm,
  // which is the only place that knows which estimator is appropriate for a given pfo.
  const float momentum[3] = {pPandoraPfo->GetMomentum().GetX(), pPandoraPfo->GetMomentum().GetY(),
                             pPandoraPfo->GetMomentum().GetZ()};
  reconstructedParticle.setMomentum(momentum);
  reconstructedParticle.setEnergy(pPandoraPfo->GetEnergy());
  reconstructedParticle.setMass(pPandoraPfo->GetMass());
  reconstructedParticle.setCharge(pPandoraPfo->GetCharge());
  reconstructedParticle.setPDG(pPandoraPfo->GetParticleId());
}
