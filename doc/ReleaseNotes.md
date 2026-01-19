# v00.02.00

* 2026-01-19 Thomas Madlener ([PR#25](https://github.com/key4hep/k4GaudiPandora/pull/25))
  - Use the properly named setter for the intrinsic cluster angle once it becomes available with EDM4hep v1.0

* 2025-12-03 Juan Miguel Carceller ([PR#24](https://github.com/key4hep/k4GaudiPandora/pull/24))
  - Add `catch(...)` to catch all and ignore the exceptions produced in DDCaloHitCreator when creating the hits in Pandora for muons.

* 2025-10-31 Juan Miguel Carceller ([PR#23](https://github.com/key4hep/k4GaudiPandora/pull/23))
  - Remove a parameter that is set to its default values and it's not present in the corresponding CLDConfig file, from which all the values are copied.

* 2025-10-31 AuroraPerego ([PR#20](https://github.com/key4hep/k4GaudiPandora/pull/20))
  - Implement a configurable Gaussian smearing for the calorimeter hits time with a fixed resolution.

* 2025-10-23 Juan Miguel Carceller ([PR#22](https://github.com/key4hep/k4GaudiPandora/pull/22))
  - Keep the property names as they were originally for `ECAL_apply_realistic_digi` and `HCAL_apply_realistic_digi`

* 2025-10-10 Thomas Madlener ([PR#21](https://github.com/key4hep/k4GaudiPandora/pull/21))
  - Make sure to run pre-commit CI on an existing nightly build
  - Fix all formatting issues

* 2025-09-26 Juan Miguel Carceller ([PR#19](https://github.com/key4hep/k4GaudiPandora/pull/19))
  - Use the `inputLocations` function in `DDCaloDigi` instead of the member since the member may be removed in a coming pull request (https://github.com/key4hep/k4FWCore/pull/345) and it's a detail of the implementation.

* 2025-07-23 jmcarcell ([PR#17](https://github.com/key4hep/k4GaudiPandora/pull/17))
  - Fix defaults in DDSimpleMuonDigi, set layers to keep to be an empty vector by default to fix a crash when the default value is used. Change the name of the output collection to be in singular `MUONOutputCollection` instead of `MUONOutputCollections`.

* 2025-06-11 jmcarcell ([PR#14](https://github.com/key4hep/k4GaudiPandora/pull/14))
  - Fix CellIDs in `DDSimpleMuonDigi`, currently set to int (like in the original processor) but we use `std::uint64_t` in EDM4hep
  - Add const where possible
  - Other minor cleanups
  - Add a test checking that the output is exactly the same as when using the Marlin wrapper. Disable running Overlay in the tests since this will (sometimes) modify calorimeter hits in place, messing up the validation.

# v00.01.00

* 2025-05-28 jmcarcell ([PR#15](https://github.com/key4hep/k4GaudiPandora/pull/15))
  - Add LANGUAGES CXX to the project call to disable checking for a C compiler

* 2025-05-13 jmcarcell ([PR#13](https://github.com/key4hep/k4GaudiPandora/pull/13))
  - Fix and clean up DDCaloDigi to give the same results as the processor. For that:
    - Use `auto` instead of `int` for the cellID in EDM4hep
    - Store `CalibrHCALOther` as a double instead of as a float
  - Add a test that runs the CLD reconstruction and checks that the result is the same with the processor and the algorithm
  - Cleanup when possible
  - Add a README listing some of the changes done compared to the original processor

* 2025-03-24 Mateusz Jakub Fila ([PR#9](https://github.com/key4hep/k4GaudiPandora/pull/9))
  - Use the properties `Input` and `Output` with `IOSvc` instead of the deprecated `input` and `output`.

* 2025-03-24 Swathi Sasikumar ([PR#4](https://github.com/key4hep/k4GaudiPandora/pull/4))
  - Port DDCaloDigi from MarlinPandora to Gaudi

* 2024-11-20 Thomas Madlener ([PR#8](https://github.com/key4hep/k4GaudiPandora/pull/8))
  - Make sure to find DD4hep first in CMake config to avoid finding conflicting versions of python.

* 2024-09-10 jmcarcell ([PR#5](https://github.com/key4hep/k4GaudiPandora/pull/5))
  - Use the Key4hepConfig flag to set the standard, compiler flags and rpath magic.

* 2024-09-09 jmcarcell ([PR#7](https://github.com/key4hep/k4GaudiPandora/pull/7))
  - Remove the old workflow for builds since we have a newer one

* 2024-08-08 Katerina Kostova ([PR#3](https://github.com/key4hep/k4GaudiPandora/pull/3))
  - Update on porting the DDSimpleMuonDigi algorithm from Marlin to Gaudi - optimised and validated 
  - Add a steering file for DDSimpleMuonDigi

* 2024-08-08 SwathiSasikumar ([PR#2](https://github.com/key4hep/k4GaudiPandora/pull/2))
  - Add initial commits with history from DDMarlinPandora and changes from Swathi

# v00-01

* This file is also automatically populated by the tagging script