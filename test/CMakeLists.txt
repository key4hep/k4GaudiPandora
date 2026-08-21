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

add_test(NAME "clone_CLDConfig"
  COMMAND bash -c "if ! test -f CLDReconstruction.py; then git clone https://github.com/key4hep/CLDConfig.git --depth 1 && cp -r CLDConfig/CLDConfig/* . && sed -i -e '/Overlay.Overlay/d' -e '/HighLevelReco.RecoMCTruthLink/d' CLDReconstruction.py; fi"
)

add_test(NAME "run_ddsim"
  COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml -G --numberOfEvents 3 --outputFile sim.edm4hep.root --gun.energy '10*GeV' --gun.particle 'mu-' --gun.distribution uniform --gun.multiplicity 10"
  )

add_test(NAME "run_ddsim_ttbar"
  COMMAND bash -c "ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml --numberOfEvents 3 --outputFile sim_ttbar.edm4hep.root --inputFile CLDConfig/test/yyxyev_000.stdhep"
  )
set_tests_properties("run_ddsim_ttbar" PROPERTIES DEPENDS "clone_CLDConfig")

add_test(NAME "run_marlin_wrapper"
  COMMAND k4run CLDReconstruction.py --inputFiles sim.edm4hep.root
)
set_tests_properties("run_marlin_wrapper" PROPERTIES DEPENDS "clone_CLDConfig;run_ddsim")

add_test(NAME "run_marlin_wrapper_ttbar"
  COMMAND k4run CLDReconstruction.py --inputFiles sim_ttbar.edm4hep.root --outputBasename output_ttbar
)
set_tests_properties("run_marlin_wrapper_ttbar" PROPERTIES DEPENDS "clone_CLDConfig;run_ddsim_ttbar")


add_test(NAME "run_Pandora"
  COMMAND k4run ${PROJECT_SOURCE_DIR}/k4GaudiPandora/options/runPandora.py
)
set_tests_properties("run_Pandora" PROPERTIES DEPENDS "run_marlin_wrapper")

add_test(NAME "run_Pandora_ttbar"
  COMMAND k4run ${PROJECT_SOURCE_DIR}/k4GaudiPandora/options/runPandora.py --IOSvc.Input output_ttbar_REC.edm4hep.root --IOSvc.Output output_pandora_ttbar.root
)
set_tests_properties("run_Pandora_ttbar" PROPERTIES DEPENDS "run_marlin_wrapper_ttbar")

add_test(NAME "compare_Pandora"
  COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/compare-pfos.py"
)
set_tests_properties("compare_Pandora" PROPERTIES DEPENDS "run_Pandora")

add_test(NAME "compare_Pandora_ttbar"
  COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/compare-pfos.py" --gaudi-file output_pandora_ttbar.root --marlin-file output_ttbar_REC.edm4hep.root
)
set_tests_properties("compare_Pandora_ttbar" PROPERTIES DEPENDS "run_Pandora_ttbar")



add_test(NAME "run_DDCaloDigi"
  COMMAND k4run ${PROJECT_SOURCE_DIR}/k4GaudiPandora/options/runDDCaloDigi.py
)
set_tests_properties("run_DDCaloDigi" PROPERTIES DEPENDS "run_marlin_wrapper")

add_test(NAME "compare_DDCaloDigi"
  COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/compare-calohits.py" --relation-collections-are-merged
)
set_tests_properties("compare_DDCaloDigi" PROPERTIES DEPENDS "run_DDCaloDigi")


add_test(NAME "run_DDSimpleMuonDigi"
  COMMAND k4run ${PROJECT_SOURCE_DIR}/k4GaudiPandora/options/runDDSimpleMuonDigi.py
)
set_tests_properties("run_DDSimpleMuonDigi" PROPERTIES DEPENDS "run_marlin_wrapper")

add_test(NAME "compare_DDSimpleMuonDigi"
  COMMAND "${PROJECT_SOURCE_DIR}/test/scripts/compare-calohits.py" --gaudi-collections GaudiMUON --marlin-collections MUON --gaudi-file output_muon_digi.root  --marlin-file output_REC.edm4hep.root --gaudi-relation GaudiRelationMUON --marlin-relation RelationMuonHit --relation-collections-are-merged
)
set_tests_properties("compare_DDSimpleMuonDigi" PROPERTIES DEPENDS "run_DDSimpleMuonDigi")

# Set the environment for all tests in this directory at once
# This follows the k4-project-template setup
get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(test_environment "\
LD_LIBRARY_PATH=\
${PROJECT_BINARY_DIR}:\
${PROJECT_BINARY_DIR}/${PROJECT_NAME}:\
${CMAKE_CURRENT_BINARY_DIR}:\
${CMAKE_CURRENT_BINARY_DIR}/genConfDir/${CURRENT_DIR_NAME}:\
$ENV{LD_LIBRARY_PATH};\
\
GAUDI_PLUGIN_PATH=\
${PROJECT_BINARY_DIR}:\
${PROJECT_BINARY_DIR}/${PROJECT_NAME}:\
${CMAKE_CURRENT_BINARY_DIR}:\
${CMAKE_CURRENT_BINARY_DIR}/genConfDir/${CURRENT_DIR_NAME}:\
$ENV{GAUDI_PLUGIN_PATH}:\
$ENV{LD_LIBRARY_PATH};\
\
PYTHONPATH=\
${PROJECT_BINARY_DIR}/${PROJECT_NAME}/genConfDir:\
${CMAKE_CURRENT_BINARY_DIR}/genConfDir:\
$ENV{PYTHONPATH};\
")
get_property(test_names DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
set_tests_properties(${test_names} PROPERTIES ENVIRONMENT "${test_environment}")
