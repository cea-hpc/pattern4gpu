set(EXAMPLE_NAME Pattern4GPU)
project(Pattern4GPU LANGUAGES C CXX)
include(../samples.common.cmake)

add_executable(Pattern4GPU Pattern4GPUModule.cc Pattern4GPUCartesian.cc main Pattern4GPU_axl.h)

# Module GeomEnv
add_library(geomenv geomenv/GeomEnvModule.cc)
target_include_directories(geomenv PUBLIC .)
target_link_libraries(geomenv PUBLIC arcane_core)

# Cartesian
add_library(cartesian cartesian/CartesianConnectivity.cc
                      cartesian/CartesianMesh.cc
                      cartesian/CellDirectionMng.cc
                      cartesian/FaceDirectionMng.cc
                      cartesian/NodeDirectionMng.cc
                      cartesian/CartesianItemSorter.cc)
target_include_directories(cartesian PUBLIC .)
target_link_libraries(cartesian PUBLIC arcane_core)

target_link_libraries(Pattern4GPU PRIVATE geomenv cartesian)

# Les axl
arcane_generate_axl(Pattern4GPU)
arcane_target_add_axl(geomenv geomenv/GeomEnv)

configure_file(Pattern4GPU.config ${CMAKE_CURRENT_BINARY_DIR} @ONLY)
#configure_file(${EXAMPLE_NAME}.arc ${CMAKE_CURRENT_BINARY_DIR} @ONLY)
arcane_add_arcane_libraries_to_target(Pattern4GPU)
target_include_directories(Pattern4GPU PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})

# Commande pour lancer via 'arcane_run'
set(RUN_COMMAND ${ARCANE_PREFIX_DIR}/bin/arcane_run -E ${CMAKE_CURRENT_BINARY_DIR}/Pattern4GPU)


# Pour les test NONREG
# Test name is also test directory
set( TEST_DIR "${CMAKE_CURRENT_LIST_DIR}" )

