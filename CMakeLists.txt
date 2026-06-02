cmake_minimum_required(VERSION 3.12)
project(DD4hepTutorials)

include(GNUInstallDirs)

find_package(DD4hep REQUIRED)

dd4hep_set_compiler_flags()
dd4hep_configure_output()

add_compile_options(-Wno-unused-variable)

add_subdirectory(detectorSegmentations)
add_subdirectory(Grainita_ECAL)
add_subdirectory(GrainitaCaloSDAction)
