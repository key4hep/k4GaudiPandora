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
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)


project(k4GaudiPandora)

# please keep this layout for version setting, used by the automatic tagging script
set(PACKAGE_VERSION_MAJOR 1)
set(PACKAGE_VERSION_MINOR 2)
set(PACKAGE_VERSION_PATCH 3)

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
FILE( GLOB _plugin_sources src/*.c* )
gaudi_add_module(k4GaudiPandoraPlugins
                SOURCES         ${_plugin_sources}
                LINK            Gaudi::GaudiKernel
                                EDM4HEP::edm4hep
				${PandoraSDK_LIBRARIES}
				${LCConent_LIBRARIES}
				DD4hep::DDCore
                                ROOT::Hist
                                ${ROOT_LIBRARIES}
                                k4FWCore::k4FWCore )

INSTALL( TARGETS k4GaudiPandoraPlugins
        EXPORT ${CMAKE_PROJECT_NAME}Targets
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/@{CMAKE_PROJECT_NAME}"
        COMPONENT dev )

TARGET_INCLUDE_DIRECTORIES(k4GaudiPandoraPlugins
                PRIVATE ./include
                        ${CMAKE_CURRENT_BINARY_DIR}
                        ${Gaudi_INCLUDE_DIRS}
			${PandoraSDK_INCLUDE_DIRS}
			${DD4hep_INCLUDE_DIRS}
			${LCContent_INCLUDE_DIRS}
                        ${ROOT_INCLUDE_DIRS})

SET( GAUDI_GENCONF_DIR "genConfDir" )
if (${Gaudi_VERSION} VERSION_LESS 35.1)
        SET( GAUDI_GENCONF_DIR "genConf" )
endif()

include(GNUInstallDirs)

include(CTest)

include(cmake/CreateProjectConfig.cmake)
