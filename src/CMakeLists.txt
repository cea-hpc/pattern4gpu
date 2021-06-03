set(EXAMPLE_NAME Pattern4GPU)
project(Pattern4GPU LANGUAGES C CXX)
include(../samples.common.cmake)

#target_sources(${EXAMPLE_NAME} PRIVATE MonFichier.cc)

# Commande pour lancer via 'arcane_run'
set(RUN_COMMAND ${ARCANE_PREFIX_DIR}/bin/arcane_run -E ${CMAKE_CURRENT_BINARY_DIR}/Pattern4GPU)


# Pour les test NONREG
# Test name is also test directory
set( TEST_DIR "${CMAKE_CURRENT_LIST_DIR}" )

