cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
set(BUILD_SHARED_LIBS TRUE)
if (NOT Arcane_DIR)
#  include(ArcaneCompilerConfig.cmake)
  message(STATUS "ArcaneRoot = ${Arcane_DIR}")
  message(STATUS "AxlstarRoot = ${Axlstar_DIR}")
  message(STATUS "ArcconRoot = ${Arccon_DIR}")
  message(STATUS "ArccoreRoot = ${Arccore_DIR}")
endif()
project(ArcaneSamples LANGUAGES CXX)

enable_testing()

add_subdirectory(src)
