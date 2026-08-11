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

# The track state at the calorimeter face can come from two places:
#
#  - OFF (default): the legacy pipeline. k4Reco::GaudiTrkUtils re-propagates the
#    track with DDKalTest to the ECal endcap face and an alternative track state
#    is appended, so that Pandora can try both it and the one already stored in
#    the EDM. This exists because GaudiTrkUtils::createTrackStateAtCaloFace does
#    not navigate to the calorimeter: it propagates to the ECal barrel and to
#    the ECal endcap and then keeps whichever intersection is closest to the
#    last tracker hit, which picks the wrong face in the barrel/endcap
#    transition region. It pulls in LCIO, KalTest and DDKalTest.
#
#  - ON: the track state stored in the EDM is used as is. k4ActsTracking's
#    CKFTrackingAlg (ExtrapolateToCalo=True) obtains it by propagating with an
#    Acts::Navigator over the calorimeter face surfaces, so the face is resolved
#    by navigation rather than guessed and the alternative state is redundant.
#
# Note that this does not add a build dependency on k4ActsTracking; it is a
# statement about the input data, not about what we link against.
option(K4GAUDIPANDORA_USE_K4ACTSTRACKING_TRACKS
       "Take the track state at the calorimeter from the input EDM instead of re-propagating it with DDKalTest. Requires input tracks produced by k4ActsTracking's CKFTrackingAlg with ExtrapolateToCalo=True, and drops the dependency on k4Reco::GaudiTrkUtils (and hence on LCIO, KalTest and DDKalTest)."
       OFF)

if(NOT K4GAUDIPANDORA_USE_K4ACTSTRACKING_TRACKS AND NOT TARGET k4Reco::GaudiTrkUtils)
  message(FATAL_ERROR
    "k4Reco::GaudiTrkUtils was not found. It is provided by k4Reco built with "
    "BUILD_TRACKING=ON.\n"
    "To build against a stack without LCIO, configure with "
    "-DK4GAUDIPANDORA_USE_K4ACTSTRACKING_TRACKS=ON.")
endif()

message(STATUS "K4GAUDIPANDORA_USE_K4ACTSTRACKING_TRACKS: ${K4GAUDIPANDORA_USE_K4ACTSTRACKING_TRACKS}")

#---------------------------------------------------------------

include(cmake/Key4hepConfig.cmake)

include(GNUInstallDirs)


add_subdirectory(k4GaudiPandora)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

include(cmake/CreateProjectConfig.cmake)
