cmake_minimum_required(VERSION 3.16)

project(externaldftbpmodel VERSION 0.0.1 LANGUAGES C DESCRIPTION
"Toy external model for the DFTB+ code")

option(BUILD_SHARED_LIBS "Whether the library should be a shared one" TRUE)

set(sources src/externalmodel.c)

add_library(ExternalModel ${sources})

set_target_properties(ExternalModel PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION 1
  PUBLIC_HEADER include/externalmodel.h)

set_target_properties(ExternalModel PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
