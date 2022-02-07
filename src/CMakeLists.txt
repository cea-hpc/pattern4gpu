set(EXAMPLE_NAME Pattern4GPU)
project(Pattern4GPU LANGUAGES C CXX)
include(../samples.common.cmake)

# Pour MPI
find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
# fin MPI

add_executable(Pattern4GPU main.cc)
target_link_libraries(Pattern4GPU PUBLIC libpattern4gpu arcane_core)

add_library(libpattern4gpu Pattern4GPUModule.cc
                           Pattern4GPUupdateTensor.cc 
                           Pattern4GPUTestCartesian.cc 
                           Pattern4GPUBenchCartesian.cc 
                           Pattern4GPUCartesian.cc 
                           Pattern4GPUEnvOrder.cc 
                           Pattern4GPUMultiEnv.cc 
                           Pattern4GPU_axl.h)
target_include_directories(libpattern4gpu PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})

# Module GeomEnv
add_library(libgeomenv geomenv/GeomEnvModule.cc)
target_include_directories(libgeomenv PUBLIC .)
target_link_libraries(libgeomenv PUBLIC arcane_core)

# Cartesian
add_library(libcartesian cartesian/CartesianConnectivity.cc
                         cartesian/CartesianMesh.cc
                         cartesian/CellDirectionMng.cc
                         cartesian/FaceDirectionMng.cc
                         cartesian/NodeDirectionMng.cc
                         cartesian/CartesianItemSorter.cc)
target_include_directories(libcartesian PUBLIC .)
target_link_libraries(libcartesian PUBLIC arcane_core)

# Services communs pour les accélérateurs
add_library(libaccenv accenv/AccEnvDefaultService.cc)
target_include_directories(libaccenv PUBLIC .)
# Pour HWLOC (en attendant un FindHWLOC.cmake)
set(ENV_HWLOC_INCDIR "$ENV{HWLOC_INCDIR}")
set(ENV_HWLOC_LIBDIR "$ENV{HWLOC_LIBDIR}")
if (WANT_CUDA AND ENV_HWLOC_INCDIR AND ENV_HWLOC_LIBDIR)
  target_compile_definitions(libaccenv PRIVATE ACCENV_HWLOC)
  target_include_directories(libaccenv PRIVATE "${ENV_HWLOC_INCDIR}")
  target_link_libraries(libaccenv PRIVATE "-L${ENV_HWLOC_LIBDIR} -lhwloc")
  target_link_libraries(libaccenv PRIVATE cuda)
endif()
# fin HWLOC
target_link_libraries(libaccenv PUBLIC arcane_core)

# Message Passing
add_library(libmsgpass msgpass/VarSyncMng.cc 
                       msgpass/GlobalSynchronize.cc
                       msgpass/GlobalSynchronizeQueue.cc
                       msgpass/GlobalSynchronizeQueueEvent.cc
                       msgpass/GlobalSynchronizeDevThr.cc
                       msgpass/GlobalSynchronizeDevQueues.cc
                       msgpass/IncompleteGlobalSynchronizeQueue.cc
                       msgpass/PackTransfer.cc
                       msgpass/SyncItems.cc
                       msgpass/SyncBuffers.cc)
target_include_directories(libmsgpass PUBLIC .)
target_link_libraries(libmsgpass PUBLIC arcane_core)

if(WANT_KOKKOS)
target_compile_options(libpattern4gpu PUBLIC -fopenmp)
endif()
target_link_libraries(libpattern4gpu PUBLIC libgeomenv libcartesian libaccenv libmsgpass)

# Pour la gestion des accélérateurs et des sources à compiler avec le
# compilateur pour accélérateur
if(WANT_CUDA)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set_source_files_properties(
      Pattern4GPUModule.cc 
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
      msgpass/VarSyncMng.cc
      msgpass/GlobalSynchronize.cc
      msgpass/GlobalSynchronizeQueue.cc
      msgpass/GlobalSynchronizeQueueEvent.cc
      msgpass/GlobalSynchronizeDevThr.cc
      msgpass/GlobalSynchronizeDevQueues.cc
      msgpass/IncompleteGlobalSynchronizeQueue.cc
      msgpass/PackTransfer.cc
      msgpass/SyncItems.cc
      msgpass/SyncBuffers.cc
      PROPERTIES LANGUAGE CUDA)
  target_link_libraries(libpattern4gpu PUBLIC arcane_accelerator_cuda_runtime)
  target_link_libraries(libgeomenv PUBLIC arcane_accelerator_cuda_runtime)
  target_link_libraries(libcartesian PUBLIC arcane_accelerator_cuda_runtime)
  target_link_libraries(libaccenv PUBLIC arcane_accelerator_cuda_runtime)
  target_link_libraries(libmsgpass PUBLIC arcane_accelerator_cuda_runtime)
endif()

# Kokkos
# a regarder pour transformer ca en external project
if(WANT_KOKKOS)
  find_package(Kokkos COMPONENTS separable_compilation REQUIRED)
  add_library(libpattern4gpu4kokkos Pattern4GPU4Kokkos.cc)
  kokkos_compilation(TARGET libpattern4gpu4kokkos)
  target_compile_definitions(libpattern4gpu4kokkos PUBLIC KOKKOS_IMPL_HALF_TYPE_DEFINED WITH_KOKKOS)
  target_compile_options(libpattern4gpu4kokkos PUBLIC --expt-relaxed-constexpr -g)  # ajouter -G pour le debug
  target_link_libraries(libpattern4gpu4kokkos PUBLIC Kokkos::kokkos arcane_core)
  target_link_libraries(libpattern4gpu PRIVATE libpattern4gpu4kokkos)
endif()


# Pour la gestion du profiling sur accélérateur
option(WANT_PROF_ACC "Active le profiling des points d'entrées pour accélérateur" FALSE)
if (WANT_PROF_ACC)
  target_compile_definitions(libpattern4gpu PRIVATE PROF_ACC)
  target_compile_definitions(libgeomenv PRIVATE PROF_ACC)
  target_compile_definitions(libcartesian PRIVATE PROF_ACC)
  target_compile_definitions(libaccenv PRIVATE PROF_ACC)
  target_compile_definitions(libmsgpass PRIVATE PROF_ACC)
endif()

# Les axl
arcane_generate_axl(Pattern4GPU)
arcane_target_add_axl(libgeomenv geomenv/GeomEnv)
arcane_target_add_axl(libaccenv accenv/AccEnvDefault)

configure_file(Pattern4GPU.config ${CMAKE_CURRENT_BINARY_DIR} @ONLY)
#configure_file(${EXAMPLE_NAME}.arc ${CMAKE_CURRENT_BINARY_DIR} @ONLY)
arcane_add_arcane_libraries_to_target(Pattern4GPU)
target_include_directories(Pattern4GPU PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})

# Commande pour lancer via 'arcane_run'
set(RUN_COMMAND ${ARCANE_PREFIX_DIR}/bin/arcane_run -E ${CMAKE_CURRENT_BINARY_DIR}/Pattern4GPU)


# Pour les test NONREG
# Test name is also test directory
set( TEST_DIR "${CMAKE_CURRENT_LIST_DIR}" )

