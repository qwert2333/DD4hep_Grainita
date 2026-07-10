#################################################################################
## Package : DigiCalo
#################################################################################

project(DigiCalo)

find_package(ROOT COMPONENTS RIO Tree Hist Graf Physics REQUIRED)
find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(DD4hep REQUIRED)
find_package(sipm REQUIRED)

gaudi_add_module(CaloDigi
  SOURCES
    src/SimpleSiPMDigiAlg.cpp
    src/SimpleCalibAlg.cpp
    src/SimulateSiPMwithEdep.cpp
  LINK
    k4FWCore::k4FWCore
    k4FWCore::k4Interface
    Gaudi::GaudiKernel
    EDM4HEP::edm4hep
    DD4hep::DDCore
    DD4hep::DDRec
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::Graf
    ROOT::Physics
    sipm::sipm
)

install(TARGETS CaloDigi
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
