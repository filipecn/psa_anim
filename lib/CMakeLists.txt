# ##############################################################################
# psa_anim lib                                                                ##
# ##############################################################################
set(PSA_ANIM_HEADERS
  ${PSA_ANIM_ROOT}/lib/psa_anim/block_mesh.h
  ${PSA_ANIM_ROOT}/lib/psa_anim/colors.h
  ${PSA_ANIM_ROOT}/lib/psa_anim/mesh_utils.h)

set(PSA_ANIM_SOURCES
  ${PSA_ANIM_ROOT}/lib/psa_anim/block_mesh.cpp
  ${PSA_ANIM_ROOT}/lib/psa_anim/colors.cpp
  ${PSA_ANIM_ROOT}/lib/psa_anim/post_process.cpp
  ${PSA_ANIM_ROOT}/lib/psa_anim/mesh_utils.cpp
  ${PSA_ANIM_ROOT}/lib/psa_anim/terrain_extruder.cpp)

add_library(psa_anim STATIC ${PSA_ANIM_SOURCES} ${PSA_ANIM_HEADERS})
add_dependencies(psa_anim hermes)
target_include_directories(
  psa_anim PRIVATE  ${EXT_INCLUDES} 
                    ${PSA_ANIM_ROOT}/lib/psa_anim
                    ${HERMES_INCLUDES})

# ##############################################################################
# python                                                                      ##
# ##############################################################################
pybind11_add_module(psa_anim_py psa_anim_py.cpp)

add_dependencies(psa_anim_py psa_anim hermes)
target_link_libraries(psa_anim_py PUBLIC 
  hermes
  psa_anim
  ${HERMES_LIBRARIES}
  Eigen3::Eigen)

target_include_directories(psa_anim_py PUBLIC 
  ${HERMES_INCLUDES}
  ${PSA_ANIM_ROOT}/lib/psa_anim)


