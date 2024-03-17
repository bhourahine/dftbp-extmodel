cmake_minimum_required(VERSION 3.16)
project(ExternalModel VERSION 0.0.1 DESCRIPTION "Toy external model for the DFTB+ code")

set(sources src/externalmodel.c)

add_library(externalmodel_objlib OBJECT ${sources})

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(externalmodel_objlib PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION 1
  PUBLIC_HEADER include/externalmodel.h)

set_target_properties(externalmodel_objlib PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

target_include_directories(externalmodel_objlib PUBLIC ${includedir})
