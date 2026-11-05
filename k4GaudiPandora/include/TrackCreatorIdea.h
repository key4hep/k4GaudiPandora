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
#ifndef TrackCreatorIdea_h
#define TrackCreatorIdea_h 1

#include "DDTrackCreatorBase.h"

/**
 *  @brief  Track creator for the IDEA detector.
 *
 *  The tracks arrive from TracksFromGenParticles already selected (calo-reaching, no ghost helices)
 *  and already extrapolated to the calorimeter face, so this creator neither applies quality cuts
 *  nor builds the DDKalTest tracking system: only CreateTracks differs from the base.
 */
class TrackCreatorIdea : public DDTrackCreatorBase {
public:
  TrackCreatorIdea(const Settings& settings, pandora::Pandora& pandora, const Gaudi::Algorithm* alg);
  ~TrackCreatorIdea() override = default;

  pandora::StatusCode CreateTracks(const std::vector<edm4hep::Track>& tracks) override;
};

#endif
