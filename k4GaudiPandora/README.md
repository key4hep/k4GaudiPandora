<!--
Copyright (c) 2020-2024 Key4hep-Project.

This file is part of Key4hep.
See https://key4hep.github.io/key4hep-doc/ for further info.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-->
## Build options

### `K4GAUDIPANDORA_USE_DDKALTEST` (default `ON`)

Controls whether tracks are extrapolated to the calorimeter face internally, with DDKalTest through
`k4Reco::GaudiTrkUtils`.

With the default `ON`, the legacy DDMarlinPandora behaviour is built: for tracks whose
`edm4hep::TrackState::AtCalorimeter` state is in the barrel,
`DDTrackCreatorBase::GetTrackStatesAtCalo` re-propagates the track to the ECal endcap face and
passes the result to pandora as an *additional* track state. LCContent's
`TrackClusterAssociationAlgorithm` then tries both. This matters because a particle entering
through the barrel can still shower in the endcap, so which face is the relevant one is not known
until the cluster is (see
[iLCSoft/DDMarlinPandora#12](https://github.com/iLCSoft/DDMarlinPandora/pull/12)). It requires
`k4Reco::GaudiTrkUtils`, and therefore LCIO, KalTest and DDKalTest.

With `OFF`, only the track states already present in the input EDM are passed on, and
`k4Reco::GaudiTrkUtils` — and with it LCIO — is no longer needed. This requires the extrapolation to
the calorimeter to have been done upstream: `k4ActsTracking`'s `CKFTrackingAlg` does it with
`ExtrapolateToCalo=True`. Note that only the states present on the track are used, so if the
upstream extrapolation stores a single one, only that face is offered to pandora.

Note that `OFF` does not make k4GaudiPandora depend on k4ActsTracking; nothing is linked or included
from it.

Caveats when `OFF`:

- The input tracks **must** carry an `AtCalorimeter` track state. Tracks without one are dropped by
  the track creators, which report `Failed to extract a track`.
- `TrackStateTolerance` has no effect, since it only bounds the acceptance radius of the endcap
  state. `TrackSystemName` has no effect either — but note it is already inert in the `ON` build
  too, as `GaudiDDKalTest` is used regardless of its value.

DDCaloDigi has been ported. The following changes have been done:
- The function `getLayerConfig` has been included inside `initialize()` to avoid
  having it run multiple times for every hit in the input. The member
  `m_layerTypes` is being used instead.
- The functions `digitalEcalCalibCoeff` and `analogueEcalCalibCoeff` have been
  merged into `ecalCalibCoeff` since they were the same function.
- Time smearing: By default in `DDCaloDigi.cc` the simhits time is taken "as is"
  when computing the rechits time. There is now the possibility to change this
  behaviour and get more realistic rechits using the `enableHitsTimeSmearing`
  flag. If set to `True`, it applies a Gaussian smearing to the simhits time.
  The `sigma` of the Gaussian is configurable as well with the
  `{E/H}CALTimeResolution` flag (in `ns`).

DDMarlinPandora has been ported
- DDPfoCreator: `SetRecoParticleReferencePoint` has been removed
