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
#include "TrackCreatorIdea.h"

#include "Pandora/PandoraEnumeratedTypes.h"
#include "Pandora/PandoraInputTypes.h"

#include <Gaudi/Algorithm.h>

#include <cmath>
#include <limits>

TrackCreatorIdea::TrackCreatorIdea(const Settings& settings, pandora::Pandora& pandora,
                                   const Gaudi::Algorithm* algorithm)
    : DDTrackCreatorBase(settings, pandora, algorithm) {
  // Deliberately no InitialiseTrackingSystem(): the track states arrive already extrapolated to the
  // calorimeter face from TracksFromGenParticles, so GetTrackStatesAtCalo is never called and the
  // DDKalTest tracking system is never needed.
}

pandora::StatusCode TrackCreatorIdea::CreateTracks(const std::vector<edm4hep::Track>& tracks) {
  // Track selection (calo-reaching, no ghost helices) is done upstream in
  // TracksFromGenParticles; create a Pandora track for every input track.
  for (const auto& pTrack : tracks) {
    // Wrap the whole per-track body: a single track with a bad (e.g. non-finite) parameter is
    // logged and skipped rather than aborting the event.
    try {
      // Note: copy-paste of the DDTrackCreatorCLIC
      // Take the first track state for the parameters
      const auto& trackState = pTrack.getTrackStates()[0];

      // Proceed to create the pandora track
      object_creation::TrackParameters trackParameters;
      trackParameters.m_d0 = trackState.D0;
      trackParameters.m_z0 = trackState.Z0;
      trackParameters.m_pParentAddress = &pTrack;

      const float signedCurvature = trackState.omega;
      trackParameters.m_particleId = 0; // no PID at this stage
      trackParameters.m_mass = 0.;      // no mass at this stage

      if (signedCurvature != 0.f)
        trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

      GetTrackStates(pTrack, trackParameters);
      // FIXME double-check TracksFromGenParticle's track state at calo
      float r_calo = trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitudeSquared();
      trackParameters.m_reachesCalorimeter = r_calo > std::numeric_limits<float>::epsilon() ? true : false;

      trackParameters.m_canFormPfo = true;
      trackParameters.m_canFormClusterlessPfo = true;

      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(m_pandora, trackParameters))
    } catch (pandora::StatusCodeException& statusCodeException) {
      m_algorithm.error() << "Failed to extract a track: " << statusCodeException.ToString() << endmsg;
      m_algorithm.debug() << " failed track : " << pTrack << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}
