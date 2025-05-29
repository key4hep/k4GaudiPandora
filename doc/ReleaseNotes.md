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