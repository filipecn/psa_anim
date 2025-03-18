# ##############################################################################
# CMAKE                                                                        #
# ##############################################################################
cmake_minimum_required(VERSION 3.12..3.30 FATAL_ERROR)

project(
  PSA_ANIM
  VERSION 0.01
  DESCRIPTION "Original code of the paper 'Digital Animation of Powder-Snow Avalanches' [SIGGRAPH 2025] "
  LANGUAGES CXX)

set(PSA_ANIM_ROOT ${CMAKE_CURRENT_LIST_DIR})

# ##############################################################################
# OPTIONS                                                                      #
# ##############################################################################
option(BUILD_PSA_ANIM_TESTS "build unit tests" OFF)

# ##############################################################################
# COMPILATION                                                                  #
# ##############################################################################
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE "ON")

if (CMAKE_COMPILER_IS_GNUCXX)
    set(DEBUG_FLAGS "${DEBUG_FLAGS} -g -pg -Wall -Wextra -O0 -fprofile-arcs -ftest-coverage --coverage -fPIC")
    set(RELEASE_FLAGS "-fopenmp -O3 -fPIC")
endif (CMAKE_COMPILER_IS_GNUCXX)

set(CMAKE_CXX_FLAGS ${RELEASE_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_CONFIGURATION_TYPES Release)

# ##############################################################################
# DEPENDENCIES                                                                 #
# ##############################################################################
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenVDB REQUIRED)

include(ExternalProject)
list(APPEND CMAKE_MODULE_PATH "${PSA_ANIM_ROOT}/ext")
include(hermes)
include(pybind11)
include(lodepng)
set(EXT_INCLUDES "${PSA_ANIM_ROOT}/ext")

# ##############################################################################
# SOURCE                                                                      ##
# ##############################################################################
add_subdirectory(lib)
set(PSA_ANIM_INCLUDES 
  ${PSA_ANIM_ROOT}/lib
  ${PSA_ANIM_ROOT}/lib/psa_anim
  ${HERMES_INCLUDES}
  ${EXT_INCLUDES})
set(PSA_ANIM_LIBRARIES psa_anim hermes)

# ##############################################################################
# TOOLS                                                                       ##
# ##############################################################################
add_subdirectory(tools)

# ##############################################################################
# psa_anim_rc                                                                 ##
# ##############################################################################
set(PSA_ANIM_PY_PATH ${CMAKE_BINARY_DIR}/lib)
set(TOOLS_PATH ${CMAKE_BINARY_DIR}/tools)
set(SCRIPTS_PATH ${PSA_ANIM_ROOT}/scripts)
set(TUTORIALS_PATH ${PSA_ANIM_ROOT}/tutorials)
set(RENDER_TOOLS_PATH ${PSA_ANIM_ROOT}/render)
configure_file(psa_anim_variables.sh.in psa_anim_variables.sh @ONLY)

# ##############################################################################
# TESTS                                                                       ##
# ##############################################################################
if(BUILD_PSA_ANIM_TESTS)
  include(catch2)
  add_subdirectory(tests)
endif(BUILD_PSA_ANIM_TESTS)
