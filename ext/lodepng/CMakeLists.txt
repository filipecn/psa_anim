project(LODEPNG)

cmake_minimum_required(VERSION 2.8.7)

add_library(lodepng lodepng.h lodepng.cpp)
set_target_properties(lodepng PROPERTIES
	OUTPUT_NAME "lodepng"
  FOLDER "LODEPNG")

