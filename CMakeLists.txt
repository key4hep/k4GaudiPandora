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
CMAKE_MINIMUM_REQUIRED(VERSION 3.19)


project(k4GaudiPandora LANGUAGES CXX)

# please keep this layout for version setting, used by the automatic tagging script
set(PACKAGE_VERSION_MAJOR 0)
set(PACKAGE_VERSION_MINOR 1)
set(PACKAGE_VERSION_PATCH 0)

set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

list(APPEND CMAKE_MODULE_PATH $ENV{PANDORAPFA}/cmakemodules)

# dependencies are declared like this
find_package(ROOT COMPONENTS RIO Tree)
find_package(DD4hep)
find_package(EDM4HEP)
find_package(k4FWCore)
find_package(PandoraSDK)
find_package(LCContent)
find_package(Gaudi)

#---------------------------------------------------------------

include(cmake/Key4hepConfig.cmake)

include(GNUInstallDirs)


add_subdirectory(k4GaudiPandora)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

include(cmake/CreateProjectConfig.cmake)
