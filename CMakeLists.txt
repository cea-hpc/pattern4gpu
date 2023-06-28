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

set(P4G_EXE "${ArcaneSamples_BINARY_DIR}/src/Pattern4GPU")

set(P4G_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/tests")
set(OUTPUT_DIR "${ArcaneSamples_BINARY_DIR}/output")


# Get test names and number of cores from JSON file
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/ctest/test_list.json" TESTLIST_JSON)

string(JSON NB_TEST LENGTH "${TESTLIST_JSON}" "write_vs_read_test_list")
math(EXPR NB_TEST "${NB_TEST} - 1")

foreach(X RANGE ${NB_TEST})
  string(JSON TEST_DIR GET "${TESTLIST_JSON}" "write_vs_read_test_list" ${X} "test-dir")

  string(JSON WRITE_TEST_NAME GET "${TESTLIST_JSON}" "write_vs_read_test_list" ${X} "write-test-name")
  string(JSON WRITE_ARCANE_PARAM GET "${TESTLIST_JSON}" "write_vs_read_test_list" ${X} "write-arcane-param")

  string(JSON READ_TEST_NAME GET "${TESTLIST_JSON}" "write_vs_read_test_list" ${X} "read-test-name")
  string(JSON READ_ARCANE_PARAM GET "${TESTLIST_JSON}" "write_vs_read_test_list" ${X} "read-arcane-param")

  string(JSON MAX_ITER GET "${TESTLIST_JSON}" "write_vs_read_test_list" ${X} "max-iter")
  string(JSON EPSILON GET "${TESTLIST_JSON}" "write_vs_read_test_list" ${X} "epsilon")

  set(TEST_NAME "${WRITE_TEST_NAME}:compare_WRITE_and_READ")

  message(STATUS "P4G -- WRITE_VS_READ: WRITE_TEST_NAME=${WRITE_TEST_NAME} READ_TEST_NAME=${READ_TEST_NAME} MAX_ITER=${MAX_ITER}")

  configure_file(ctest/launch_test.sh.in ${CMAKE_CURRENT_BINARY_DIR}/test_scripts/write_vs_read/launch_${WRITE_TEST_NAME}.sh @ONLY)
  add_test(NAME ${TEST_NAME} COMMAND /bin/sh ${CMAKE_CURRENT_BINARY_DIR}/test_scripts/write_vs_read/launch_${WRITE_TEST_NAME}.sh)

endforeach()
