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

Selects where the track states at the calorimeter faces that pandora needs are computed.
The default `ON` preserves the legacy code behaviour. Setting this to `OFF` removes the
transitive dependencies to LCIO, KalTest and DDKalTest.

With the default `ON`, they are **recomputed here with DDKalTest**. For tracks whose
`edm4hep::TrackState::AtCalorimeter` state is in the barrel,
`DDTrackCreatorBase::GetTrackStatesAtCalo` recomputes the extrapolation to the ECal endcap face and
passes it to pandora as an *additional* track state, so that LCContent's
`TrackClusterAssociationAlgorithm` can consider both. A particle entering through the barrel can
still shower in the endcap, so which of the two faces is the relevant one is not known until the
cluster is (see
[iLCSoft/DDMarlinPandora#12](https://github.com/iLCSoft/DDMarlinPandora/pull/12)). This uses
`k4Reco::GaudiTrkUtils` (the source of the dependency on LCIO, KalTest and DDKalTest).

With `OFF`, they are expected to have been **computed upstream**, as `k4ActsTracking`'s
algorithms does with `ExtrapolateToCalo=True`. Nothing is recomputed here: *all*
`AtCalorimeter` track states found on the input track are forwarded to pandora, so the number of
faces it can choose between is whatever the upstream extrapolation stored.

Note that `OFF` does not make k4GaudiPandora depend on k4ActsTracking; nothing is linked or included
from it.

Caveats when `OFF`:

- The input tracks **must** carry an `AtCalorimeter` track state. Tracks without one are dropped by
  the track creators, which report `Failed to extract a track`.
- `TrackStateTolerance` has no effect, since it only bounds the acceptance radius of the endcap
  state.

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

## Theta-energy calibration

`DDPandoraPFANewAlgorithm` can correct the energy of each cluster using a table of correction factors. The factor
depends on the cluster's polar angle (theta) and its energy. Two separate tables can be given: one for
electromagnetic (EM) clusters and one for hadronic clusters. Each table is stored in a JSON file.

### Settings

| Setting | Description |
|---|---|
| `ElectromagneticThetaEnergyCorrectionFile` | Path to the JSON file with the EM correction table. |
| `HadronicThetaEnergyCorrectionFile` | Path to the JSON file with the hadronic correction table. |

Both settings are empty by default. In that case the matching plugin (see below) is still available to Pandora, but
it does not change any energies.

A relative path is taken relative to the directory the job is run from, in the same way as `PandoraSettingsXmlFile`.

The files are read when the job starts. If a file cannot be opened, or the table in it is not valid, the job stops
and prints an error explaining the problem.

### Pandora settings XML

The corrections are applied by two Pandora plugins:

| Plugin name | Corrects | Uses the table from |
|---|---|---|
| `PhotonEMNonLinearity` | EM clusters | `ElectromagneticThetaEnergyCorrectionFile` |
| `HadronicThetaEnergyBinned` | hadronic clusters | `HadronicThetaEnergyCorrectionFile` |

A correction only takes effect if the Pandora settings XML lists its plugin:

```xml
<ElectromagneticEnergyCorrectionPlugins>PhotonEMNonLinearity</ElectromagneticEnergyCorrectionPlugins>
<HadronicEnergyCorrectionPlugins>HadronicThetaEnergyBinned</HadronicEnergyCorrectionPlugins>
```

If these entries already list other plugins, add the new names to the existing list.

The `ConeBasedMerging` and `ProximityBasedMerging` algorithms can also use the hadronic table when they compare
cluster energies with track momenta. To enable this, add the following to the settings of each of these algorithms
in the XML:

```xml
<UseThetaEnergyCorrectionForTrackComparison>true</UseThetaEnergyCorrectionForTrackComparison>
<ThetaEnergyCorrectionName>HadronicThetaEnergyBinned</ThetaEnergyCorrectionName>
```

### JSON file format

Example of an EM table with two theta bins and two energy bins:

```json
{
  "theta_edges": [0.0, 1.5708, 3.1416],
  "energy_edges": [0.0, 50.0, 1000.0],
  "scales": [1.02, 1.01, 1.03, 1.02],
  "metadata": {
    "energy_basis": "em"
  }
}
```

| Key | Description |
|---|---|
| `theta_edges` | Bin edges in theta, in radians. Theta always lies between 0 and π. |
| `energy_edges` | Bin edges in energy, in GeV. |
| `scales` | The correction factors, one per bin, as a single list. |
| `metadata.energy_basis` | `"em"` for the EM file, `"hadronic"` for the hadronic file. |

The table must meet these requirements:

- `theta_edges` and `energy_edges` each have at least two values, and each value is larger than the one before it.
- `scales` has one value for every combination of a theta bin and an energy bin, which is
  (number of theta edges − 1) × (number of energy edges − 1) values.
- The values in `scales` are grouped by theta bin: first all energy bins of the first theta bin, then all energy bins
  of the second theta bin, and so on. In the example above, the first two values belong to the first theta bin.

Any other keys in the file, such as `domain` or `counts` written by the calibration scripts, are ignored.

### Clusters outside the table

The corrected energy is the original energy multiplied by the factor for the cluster's bin. A cluster that falls
outside the table keeps its original energy, with one exception for high energies:

- Theta below the first edge, or equal to or above the last edge: not corrected.
- Energy below the first edge: not corrected.
- Energy equal to or above the last edge: corrected with the factor of the last energy bin.
