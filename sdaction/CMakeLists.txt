project(simplecalo)

file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.hh
)

add_library(sdaction SHARED ${sources} ${headers})

target_include_directories(sdaction PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(sdaction PROPERTIES PUBLIC_HEADER "${headers}")

target_link_libraries(
  sdaction
  DD4hep::DDCore
  DD4hep::DDG4
)

dd4hep_generate_rootmap(sdaction)

install(TARGETS sdaction EXPORT DetectorTargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)
