CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
project(DD4hepTutorials)

include(GNUInstallDirs)

find_package(DD4hep REQUIRED)

dd4hep_set_compiler_flags()
dd4hep_configure_output()

# Collect all source files
set(SOURCES
  ${PROJECT_SOURCE_DIR}/simplecalo1/src/simplecalo1.cpp
  ${PROJECT_SOURCE_DIR}/simplecalo2/src/simplecalo2.cpp
)
add_compile_options(-Wno-unused-variable)

# Create a single plugin
dd4hep_add_plugin(${PROJECT_NAME} SOURCES ${SOURCES} USES
  DD4hep::DDCore
  DD4hep::DDCond
  ROOT::Core
  ROOT::Geom
  ROOT::GenVector
  ROOT::MathCore
)

# Add both include directories
target_include_directories(${PROJECT_NAME} PRIVATE
  ${PROJECT_SOURCE_DIR}/simplecalo1/include
  ${PROJECT_SOURCE_DIR}/simplecalo2/include
)

add_subdirectory(sdaction)
