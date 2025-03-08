#---------------------------------------------------------------
FILE( GLOB _plugin_sources src/*.c* )
gaudi_add_module(k4GaudiPandoraPlugins
                SOURCES         ${_plugin_sources}
		                ../Utility/MarlinUtil/01-08/source/ClusterShapes.cc
                                ../Utility/MarlinUtil/01-08/source/HelixClass.cc
                                ../Utility/MarlinUtil/01-08/source/LineClass.cc
                LINK            Gaudi::GaudiKernel
                                EDM4HEP::edm4hep
				${PandoraSDK_LIBRARIES}
				${LCContent_LIBRARIES}
				${GSL_LIBRARIES}
                                ${LCIO_LIBRARIES}
                                ${MarlinTrk_LIBRARIES}
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
		PUBLIC   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>/Utility/MarlinUtil/01-08/source
                PRIVATE ./include
                        ${CMAKE_CURRENT_BINARY_DIR}
                        ${Gaudi_INCLUDE_DIRS}
			${GSL_INCLUDE_DIRS}
			${PandoraSDK_INCLUDE_DIRS}
                        ${LCConent_INCLUDE_DIRS}
                        ${LCIO_INCLUDE_DIRS}
                        ${MarlinTrk_INCLUDE_DIRS}
			${DD4hep_INCLUDE_DIRS}
			${LCContent_INCLUDE_DIRS}
                        ${ROOT_INCLUDE_DIRS})
