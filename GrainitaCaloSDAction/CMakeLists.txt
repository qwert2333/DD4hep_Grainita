project(GrainitaCaloSDAction)

file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/GrainitaCaloSDAction.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/src/*.h
)

add_library(GrainitaCaloSDAction SHARED ${sources} ${headers})

target_include_directories(GrainitaCaloSDAction PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../detectorSegmentations/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(GrainitaCaloSDAction PROPERTIES PUBLIC_HEADER "${headers}")

target_link_libraries(
  GrainitaCaloSDAction
  DD4hep::DDCore
  DD4hep::DDG4
  detectorSegmentations
)

dd4hep_generate_rootmap(GrainitaCaloSDAction)

install(TARGETS GrainitaCaloSDAction EXPORT DetectorTargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)
