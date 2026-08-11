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
set(PACKAGE_VERSION_MINOR 2)
set(PACKAGE_VERSION_PATCH 0)

set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

list(APPEND CMAKE_MODULE_PATH $ENV{PANDORAPFA}/cmakemodules)

# dependencies are declared like this
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

# Whether to extrapolate to the calorimeter internally, using DDKalTest through
# k4Reco::GaudiTrkUtils:
#
#  - ON (default): for tracks whose track state at the calorimeter is in the
#    barrel, the track is re-propagated to the ECal endcap face and the result is
#    appended as an additional track state, so that pandora can try both. This is
#    the legacy DDMarlinPandora behaviour, and it pulls in LCIO, KalTest and
#    DDKalTest.
#
#  - OFF: only the track states already present in the input EDM are used, and
#    k4Reco::GaudiTrkUtils is not needed. Whether that is enough depends on what
#    produced the tracks: it requires the extrapolation to the calorimeter to have
#    been done upstream, as k4ActsTracking's CKFTrackingAlg does with
#    ExtrapolateToCalo=True.
option(K4GAUDIPANDORA_USE_DDKALTEST
       "Extrapolate tracks to the calorimeter face internally with DDKalTest, adding a track state at the ECal endcap on top of the ones present in the input EDM. Requires k4Reco::GaudiTrkUtils, and hence LCIO, KalTest and DDKalTest."
       ON)

if(K4GAUDIPANDORA_USE_DDKALTEST AND NOT TARGET k4Reco::GaudiTrkUtils)
  message(FATAL_ERROR
    "k4Reco::GaudiTrkUtils was not found. It is provided by k4Reco built with "
    "BUILD_TRACKING=ON.\n"
    "To build against a stack without LCIO, configure with "
    "-DK4GAUDIPANDORA_USE_DDKALTEST=OFF. The track states are then taken from the "
    "input EDM as they are.")
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
