/// Copyright (c) 2020, FilipeCN.
///
/// The MIT License (MIT)
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.
///
/// \file main.cpp
/// \author FilipeCN (filipedecn@gmail.com)
/// \date 2020-02-14
///
/// \brief

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <terrain_extruder.h>

using namespace nadare;


TEST_CASE("plane") {
  hermes::Path assets_path(ASSETS_PATH);
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  auto asset_path = assets_path + "plane.stl";
  TerrainExtruder extruder;
  REQUIRE(extruder.load(asset_path));

  extruder.save("test.stl");
  REQUIRE(extruder["top"].vertices.size() == 4);
  REQUIRE(extruder["terrain"].vertices.size() == 4);
  REQUIRE(extruder["walls"].vertices.size() == 8);
  REQUIRE(extruder["inlet"].vertices.size() == 4);
  REQUIRE(extruder["outlet"].vertices.size() == 4);

  REQUIRE(extruder["top"].faces.size() == 2);
  REQUIRE(extruder["terrain"].faces.size() == 2);
  REQUIRE(extruder["walls"].faces.size() == 4);
  REQUIRE(extruder["inlet"].faces.size() == 2);
  REQUIRE(extruder["outlet"].faces.size() == 2);

}

TEST_CASE("plane45") {
  hermes::Path assets_path(ASSETS_PATH);
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  auto asset_path = assets_path + "plane45.stl";
  TerrainExtruder extruder;
  extruder.setSlopeDirection({0, 0}, {std::sqrt(2.f), std::sqrt(2.f)});
  REQUIRE(extruder.load(asset_path));

  REQUIRE(extruder["top"].vertices.size() == 4);
  REQUIRE(extruder["terrain"].vertices.size() == 4);
  REQUIRE(extruder["walls"].vertices.size() == 8);
  REQUIRE(extruder["inlet"].vertices.size() == 4);
  REQUIRE(extruder["outlet"].vertices.size() == 4);

  REQUIRE(extruder["top"].faces.size() == 2);
  REQUIRE(extruder["terrain"].faces.size() == 2);
  REQUIRE(extruder["walls"].faces.size() == 4);
  REQUIRE(extruder["inlet"].faces.size() == 2);
  REQUIRE(extruder["outlet"].faces.size() == 2);

  extruder.save("test45.stl");
}

TEST_CASE("wolfsgrube") {
  hermes::Path assets_path(ASSETS_PATH);
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  auto asset_path = assets_path + "wolfsgrube.stl";
  TerrainExtruder extruder;
  REQUIRE(extruder.load(asset_path));
  extruder.save("");
}

TEST_CASE("load patch") {
  hermes::Path assets_path(ASSETS_PATH);
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  auto asset_path = assets_path + "test.stl";
  TerrainExtruder extruder;
  REQUIRE(extruder.load(asset_path, "terrain"));
}