#[[
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
]]


find_package(k4geo REQUIRED)

function(set_test_env testname)
  set_property(TEST ${testname} APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}:${PROJECT_BINARY_DIR}/${PROJECT_NAME}:${PROJECT_BINARY_DIR}/${PROJECT_NAME}/genConfDir/${PROJECT_NAME}:${PROJECT_BINARY_DIR}/test/k4FWCoreTest:${PROJECT_BINARY_DIR}/test/k4FWCoreTest/genConfDir/k4FWCoreTest:$<$<TARGET_EXISTS:ROOT::Core>:$<TARGET_FILE_DIR:ROOT::Core>>:$<$<TARGET_EXISTS:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:EDM4HEP::edm4hep>>:$<$<TARGET_EXISTS:podio::podio>:$<TARGET_FILE_DIR:podio::podio>>:$ENV{LD_LIBRARY_PATH}")
  set_property(TEST ${testname} APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/${PROJECT_NAME}/genConfDir:$ENV{PYTHONPATH}")
endfunction()

add_test(NAME "Clone CLDConfig"
  COMMAND bash -c "git clone https://github.com/key4hep/CLDConfig.git --depth 1 && cp -r CLDConfig/CLDConfig/* ."
)
set_test_env("Clone CLDConfig")

add_test(NAME "Run ddsim"
  COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml -G --numberOfEvents 3 --outputFile sim.edm4hep.root --gun.energy '10*GeV' --gun.particle 'mu-' --gun.distribution uniform --gun.multiplicity 10"
  )
set_test_env("Run ddsim")

add_test(NAME "Run the Marlin wrapper"
  COMMAND k4run CLDConfig/CLDConfig/CLDReconstruction.py --inputFiles sim.edm4hep.root
)
set_test_env("Run the Marlin wrapper")
set_tests_properties("Run the Marlin wrapper" PROPERTIES DEPENDS "Clone CLDConfig;Run ddsim")

add_test(NAME "Run the DDCaloDigi algorithm"
  COMMAND k4run ${PROJECT_SOURCE_DIR}/k4GaudiPandora/options/runDDCaloDigi.py
)
set_test_env("Run the DDCaloDigi algorithm")
set_tests_properties("Run the DDCaloDigi algorithm" PROPERTIES DEPENDS "Run the Marlin wrapper")

add_test(NAME "Compare calorimeter hits between Marlin and Gaudi"
  COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/compare-calohits.py"
)
set_test_env("Compare calorimeter hits between Marlin and Gaudi")
set_tests_properties("Compare calorimeter hits between Marlin and Gaudi" PROPERTIES DEPENDS "Run the DDCaloDigi algorithm")

