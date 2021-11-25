set(EXAMPLE_NAME Pattern4GPU)
project(Pattern4GPU LANGUAGES C CXX)
include(../samples.common.cmake)

add_executable(Pattern4GPU Pattern4GPUModule.cc
                           Pattern4GPUupdateTensor.cc 
                           Pattern4GPUTestCartesian.cc 
                           Pattern4GPUBenchCartesian.cc 
                           Pattern4GPUCartesian.cc 
                           Pattern4GPUEnvOrder.cc 
                           Pattern4GPUMultiEnv.cc 
                           main 
                           Pattern4GPU_axl.h)

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

# Services communs pour les accélérateurs
add_library(accenv accenv/AccEnvDefaultService.cc)
target_include_directories(accenv PUBLIC .)
target_link_libraries(accenv PUBLIC arcane_core)

target_link_libraries(Pattern4GPU PRIVATE geomenv cartesian accenv)

# Pour la gestion des accélérateurs et des sources à compiler avec le
# compilateur pour accélérateur
if(WANT_CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set_source_files_properties(Pattern4GPUModule.cc 
      Pattern4GPUupdateTensor.cc 
      Pattern4GPUTestCartesian.cc 
      Pattern4GPUBenchCartesian.cc 
      Pattern4GPUCartesian.cc 
      Pattern4GPUEnvOrder.cc 
      Pattern4GPUMultiEnv.cc 
      geomenv/GeomEnvModule.cc
      cartesian/CartesianConnectivity.cc
      cartesian/CartesianMesh.cc
      cartesian/CellDirectionMng.cc
      cartesian/FaceDirectionMng.cc
      cartesian/NodeDirectionMng.cc
      cartesian/CartesianItemSorter.cc
      accenv/AccEnvDefaultService.cc
      PROPERTIES LANGUAGE CUDA)
  target_link_libraries(Pattern4GPU PUBLIC arcane_accelerator_cuda_runtime)
  target_link_libraries(geomenv PUBLIC arcane_accelerator_cuda_runtime)
  target_link_libraries(cartesian PUBLIC arcane_accelerator_cuda_runtime)
  target_link_libraries(accenv PUBLIC arcane_accelerator_cuda_runtime)
endif()

# Pour la gestion du profiling sur accélérateur
option(WANT_PROF_ACC "Active le profiling des points d'entrées pour accélérateur" FALSE)
if (WANT_PROF_ACC)
  target_compile_definitions(Pattern4GPU PRIVATE PROF_ACC)
  target_compile_definitions(geomenv PRIVATE PROF_ACC)
  target_compile_definitions(cartesian PRIVATE PROF_ACC)
  target_compile_definitions(accenv PRIVATE PROF_ACC)
endif()

# Les axl
arcane_generate_axl(Pattern4GPU)
arcane_target_add_axl(geomenv geomenv/GeomEnv)
arcane_target_add_axl(accenv accenv/AccEnvDefault)

configure_file(Pattern4GPU.config ${CMAKE_CURRENT_BINARY_DIR} @ONLY)
#configure_file(${EXAMPLE_NAME}.arc ${CMAKE_CURRENT_BINARY_DIR} @ONLY)
arcane_add_arcane_libraries_to_target(Pattern4GPU)
target_include_directories(Pattern4GPU PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})

# Commande pour lancer via 'arcane_run'
set(RUN_COMMAND ${ARCANE_PREFIX_DIR}/bin/arcane_run -E ${CMAKE_CURRENT_BINARY_DIR}/Pattern4GPU)


# Pour les test NONREG
# Test name is also test directory
set( TEST_DIR "${CMAKE_CURRENT_LIST_DIR}" )

