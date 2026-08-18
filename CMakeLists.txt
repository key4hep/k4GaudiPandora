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
set(PACKAGE_VERSION_MINOR 3)
set(PACKAGE_VERSION_PATCH 0)

set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

list(APPEND CMAKE_MODULE_PATH $ENV{PANDORAPFA}/cmakemodules)

find_package(ROOT COMPONENTS RIO Tree)
find_package(DD4hep)
find_package(GSL REQUIRED)
find_package(EDM4HEP)
find_package(k4FWCore)
find_package(PandoraSDK REQUIRED)
find_package(LCContent REQUIRED)
find_package(Gaudi)
find_package(k4Reco)

#---------------------------------------------------------------
option(K4GAUDIPANDORA_USE_DDKALTEST
       "Recompute the extrapolation to the ECal endcap face with DDKalTest, to give pandora an additional track state at the calorimeter. Requires k4Reco::GaudiTrkUtils, and hence LCIO, KalTest and DDKalTest. When OFF, the track states at the calorimeter are expected to be already present on the input tracks, and all of them are passed on."
       ON)

if(K4GAUDIPANDORA_USE_DDKALTEST AND NOT TARGET k4Reco::GaudiTrkUtils)
  message(FATAL_ERROR
    "k4Reco::GaudiTrkUtils was not found. It is provided by k4Reco built with "
    "BUILD_TRACKING=ON.\n"
    "To build against a stack without LCIO, configure with "
    "-DK4GAUDIPANDORA_USE_DDKALTEST=OFF.")
endif()

message(STATUS "K4GAUDIPANDORA_USE_DDKALTEST: ${K4GAUDIPANDORA_USE_DDKALTEST}")

#---------------------------------------------------------------

include(cmake/Key4hepConfig.cmake)

include(GNUInstallDirs)


add_subdirectory(k4GaudiPandora)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

include(cmake/CreateProjectConfig.cmake)
