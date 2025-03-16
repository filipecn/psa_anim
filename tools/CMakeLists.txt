# ubuntu: set(OPENVDB_LIBRARIES /usr/local/lib/libopenvdb.a tbb blosc pthread
# boost_iostreams /usr/lib/x86_64-linux-gnu/libz.so )

#list(APPEND CMAKE_MODULE_PATH "/usr/local/lib/cmake/OpenVDB")
#find_package(OpenVDB REQUIRED)

# ##############################################################################
# tools                                                                       ##
# ##############################################################################
function(buildTool TOOL_NAME)
  message(STATUS "Generating project file for tool ${TOOL_NAME}")
  # Main
  set(MAIN_CPP ${CMAKE_CURRENT_LIST_DIR}/${TOOL_NAME}.cpp)
  add_executable(${TOOL_NAME} ${MAIN_CPP} ${SOURCE} ${SHADERS}
                              ${TOOL_BASE_SOURCE})
  add_dependencies(${TOOL_NAME} hermes psa_anim)
  target_link_libraries(
    ${TOOL_NAME} 
    ${PSA_ANIM_LIBRARIES} 
    ${OPENVDB_LIBRARIES}
    ${LODEPNG_LIBRARIES}
    OpenVDB::openvdb 
    Eigen3::Eigen)
  target_include_directories(${TOOL_NAME} PUBLIC 
    ${LODEPNG_INCLUDES}
    ${PSA_ANIM_INCLUDES})
endfunction(buildTool)

set(TOOLS
    stl2foam
    vdb_bash
    vdb_surface
    vdb_edit
    foam2vdb
    foam2obj
    stats)

foreach(TOOL ${TOOLS})
  buildtool(${TOOL})
endforeach(TOOL)


