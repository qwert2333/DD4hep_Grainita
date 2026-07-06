#################################################################################
## Package : Grainita_ECAL
#################################################################################

project(Grainita_ECAL)

set(GRAINITA_ECAL_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Grainita_ECAL_Barrel_v01.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Grainita_ECAL_Barrel_v02.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Grainita_ECAL_Matrix.cpp
#  ${CMAKE_CURRENT_SOURCE_DIR}/src/Grainita_ECAL_Barrel_v03.cpp
)

dd4hep_add_plugin(${PROJECT_NAME} SOURCES ${GRAINITA_ECAL_SOURCES} USES
  DD4hep::DDCore
  DD4hep::DDCond
  ROOT::Core
  ROOT::Geom
  ROOT::GenVector
  ROOT::MathCore
)

target_include_directories(${PROJECT_NAME} PRIVATE
  ${CMAKE_SOURCE_DIR}/detectorSegmentations/include
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  detectorSegmentations
)
