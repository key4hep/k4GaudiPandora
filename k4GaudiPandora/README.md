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

### `K4GAUDIPANDORA_USE_K4ACTSTRACKING_TRACKS` (default `OFF`)

Selects where the track state at the calorimeter face passed to pandora comes from.

With the default `OFF`, the legacy pipeline is built: `DDTrackCreatorBase::GetTrackStatesAtCalo`
re-propagates the track to the ECal endcap face with DDKalTest and passes the result to pandora as
an *alternative* to the `edm4hep::TrackState::AtCalorimeter` state already stored in the EDM.
LCContent's `TrackClusterAssociationAlgorithm` then tries both and keeps whichever finds a cluster.

This exists because `GaudiTrkUtils::createTrackStateAtCaloFace`, which produced that stored state,
does not navigate to the calorimeter: it propagates to the ECal barrel and to the ECal endcap and
keeps whichever intersection is closest to the last tracker hit. That is not a first-intersection
criterion and it picks the wrong face in the barrel/endcap transition region, which is what the
alternative state hedges against (see
[iLCSoft/DDMarlinPandora#12](https://github.com/iLCSoft/DDMarlinPandora/pull/12)). It requires
`k4Reco::GaudiTrkUtils`, and therefore LCIO, KalTest and DDKalTest.

With `ON`, the stored state is used as is and the DDKalTest re-propagation is not built, so
`k4Reco::GaudiTrkUtils` — and with it LCIO — is no longer needed. This is correct when the input
tracks come from `k4ActsTracking`'s `CKFTrackingAlg` with `ExtrapolateToCalo=True`, which obtains
the state by propagating with an `Acts::Navigator` over the calorimeter face surfaces: the face is
resolved by navigation rather than guessed, so the alternative state is redundant.

Note that setting this does not make k4GaudiPandora depend on k4ActsTracking; it is a statement
about the input data.

Caveats when `ON`:

- The input tracks **must** carry an `AtCalorimeter` track state. Tracks without one are dropped by
  the track creators, which report `Failed to extract a track`.
- `TrackStateTolerance` has no effect, since it only bounds the acceptance radius of the endcap
  state. `TrackSystemName` has no effect either — but note it is already inert in the `OFF` build
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
