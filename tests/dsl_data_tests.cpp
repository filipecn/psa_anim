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

#include <psa_anim/openfoam_poly_mesh.h>

TEST_CASE("cube") {
  hermes::Path assets_path(ASSETS_PATH);
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  auto asset_path = assets_path + "cube/constant/polyMesh";
  PolyMesh pm;
  pm.load(asset_path);
  REQUIRE(pm.vertices.size() == 8);
  REQUIRE(pm.vertices[0] == hermes::point3(0, -0.5, 0));
  REQUIRE(pm.vertices[1] == hermes::point3(0, 0.5, 0));
  REQUIRE(pm.vertices[2] == hermes::point3(0, -0.5, 1));
  REQUIRE(pm.vertices[3] == hermes::point3(0, 0.5, 1));
  REQUIRE(pm.vertices[4] == hermes::point3(1, -0.5, 0));
  REQUIRE(pm.vertices[5] == hermes::point3(1, 0.5, 0));
  REQUIRE(pm.vertices[6] == hermes::point3(1, -0.5, 1));
  REQUIRE(pm.vertices[7] == hermes::point3(1, 0.5, 1));
  for (const auto &face : pm.faces)
    REQUIRE(face.size() == 4);
  size_t faces[6][4] = {{1, 3, 7, 5}, {0, 4, 6, 2}, {0, 1, 5, 4},
                        {2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}};
  for (size_t i = 0; i < pm.faces.size(); ++i)
    for (size_t j = 0; j < 4; ++j)
      REQUIRE(faces[i][j] == pm.faces[i][j]);
  REQUIRE(pm.owner.size() == 6);
  for (size_t i : pm.owner)
    REQUIRE(i == 0);
  REQUIRE(pm.neighbour.size() == 0);
}

TEST_CASE("points dict") {
  hermes::Path assets_path(ASSETS_PATH);
  SECTION("points1") {
    hermes::Log::addOptions(hermes::logging_options::abbreviate);
    auto asset_path = assets_path + "points1";
    SECTION("dict") {
      OpenFoamDict dict(asset_path);
      HERMES_LOG_VARIABLE(dict.nodes().size());
      for (const auto &node : dict.nodes()) {
        HERMES_LOG_VARIABLE(node.first);
      }
      REQUIRE(dict.nodes().count("FoamFile"));
      REQUIRE(dict.nodes().count("list"));
    } //
    SECTION("mesh") {
      PolyMesh pm;
      pm.readPoints(asset_path);
      REQUIRE(pm.vertices.size() == 8);
      std::vector<hermes::point3> points = {
          {0, -0.5, 0}, {0, 0.5, 0}, {0, -0.5, 1}, {0, 0.5, 1},
          {1, -0.5, 0}, {1, 0.5, 0}, {1, -0.5, 1}, {1, 0.5, 1}};
      for (size_t i = 0; i < points.size(); ++i)
        REQUIRE(points[i] == pm.vertices[i]);
    } //
  } //
  SECTION("points2") {
    hermes::Log::addOptions(hermes::logging_options::abbreviate);
    auto asset_path = assets_path + "points2";
    SECTION("dict") {
      OpenFoamDict dict(asset_path);
      HERMES_LOG_VARIABLE(dict.nodes().size());
      for (const auto &node : dict.nodes()) {
        HERMES_LOG_VARIABLE(node.first);
      }
      REQUIRE(dict.nodes().count("FoamFile"));
      REQUIRE(dict.nodes().count("list"));
    } //
    SECTION("mesh") {
      PolyMesh pm;
      pm.readPoints(asset_path);
      REQUIRE(pm.vertices.size() == 8);
      std::vector<hermes::point3> points = {
          {0, -0.5, 0}, {0, 0.5, 0}, {0, -0.5, 1}, {0, 0.5, 1},
          {1, -0.5, 0}, {1, 0.5, 0}, {1, -0.5, 1}, {1, 0.5, 1}};
      for (size_t i = 0; i < points.size(); ++i)
        REQUIRE(points[i] == pm.vertices[i]);
    } //
  } //
}

TEST_CASE("faces dict") {
  hermes::Path assets_path(ASSETS_PATH);
  SECTION("faces1") {
    hermes::Log::addOptions(hermes::logging_options::abbreviate);
    auto asset_path = assets_path + "faces1";
    SECTION("dict") {
      OpenFoamDict dict(asset_path);
      HERMES_LOG_VARIABLE(dict.nodes().size());
      for (const auto &node : dict.nodes()) {
        HERMES_LOG_VARIABLE(node.first);
      }
      REQUIRE(dict.nodes().count("FoamFile"));
      REQUIRE(dict.nodes().count("list"));
    } //
    SECTION("mesh") {
      PolyMesh pm;
      pm.readFaces(asset_path);
      REQUIRE(pm.faces.size() == 6);
      std::vector<std::vector<size_t>> faces = {{1, 3, 7, 5}, {0, 4, 6, 2},
                                                {0, 1, 5, 4}, {2, 6, 7, 3},
                                                {0, 2, 3, 1}, {4, 5, 7, 6}};
      for (size_t i = 0; i < faces.size(); ++i) {
        REQUIRE(faces[i].size() == pm.faces[i].size());
        for (size_t j = 0; j < faces[i].size(); ++j)
          REQUIRE(faces[i][j] == pm.faces[i][j]);
      }
    } //
  } //
  SECTION("faces2") {
    hermes::Log::addOptions(hermes::logging_options::abbreviate);
    auto asset_path = assets_path + "faces2";
    SECTION("dict") {
      OpenFoamDict dict(asset_path);
      HERMES_LOG_VARIABLE(dict.nodes().size());
      for (const auto &node : dict.nodes()) {
        HERMES_LOG_VARIABLE(node.first);
      }
      REQUIRE(dict.nodes().count("FoamFile"));
      REQUIRE(dict.nodes().count("list"));
    } //
    SECTION("mesh") {
      PolyMesh pm;
      pm.readFaces(asset_path);
      REQUIRE(pm.faces.size() == 6);
      std::vector<std::vector<size_t>> faces = {{1, 3, 7, 5}, {0, 4, 6, 2},
                                                {0, 1, 5, 4}, {2, 6, 7, 3},
                                                {0, 2, 3, 1}, {4, 5, 7, 6}};
      for (size_t i = 0; i < faces.size(); ++i) {
        REQUIRE(faces[i].size() == pm.faces[i].size());
        for (size_t j = 0; j < faces[i].size(); ++j)
          REQUIRE(faces[i][j] == pm.faces[i][j]);
      }
    } //
  } //
}

TEST_CASE("boundary dict") {
  hermes::Path assets_path(ASSETS_PATH);
  SECTION("boundary1") {
    hermes::Log::addOptions(hermes::logging_options::abbreviate);
    auto asset_path = assets_path + "boundary1";
    SECTION("dict") {
      OpenFoamDict dict(asset_path);
      HERMES_LOG_VARIABLE(dict.nodes().size());
      for (const auto &node : dict.nodes()) {
        HERMES_LOG_VARIABLE(node.first);
      }
      REQUIRE(dict.nodes().count("FoamFile"));
      REQUIRE(dict.nodes().count("list"));
    } //
    SECTION("mesh") {
      PolyMesh pm;
      pm.readBoundary(asset_path);
      REQUIRE(pm.patches.size() == 5);
      {
        std::string patch_name = "frontAndBack";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 2);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 0);
      } //
      {
        std::string patch_name = "terrain";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 1);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 2);
      } //
      {
        std::string patch_name = "top";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 1);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 3);
      } //
      {
        std::string patch_name = "inlet";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 1);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 4);
      } //
      {
        std::string patch_name = "outlet";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 1);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 5);
      } //
    } //
  } //
  SECTION("boundary2") {
    hermes::Log::addOptions(hermes::logging_options::abbreviate);
    auto asset_path = assets_path + "boundary2";
    SECTION("dict") {
      OpenFoamDict dict(asset_path);
      HERMES_LOG_VARIABLE(dict.nodes().size());
      for (const auto &node : dict.nodes()) {
        HERMES_LOG_VARIABLE(node.first);
      }
      REQUIRE(dict.nodes().count("FoamFile"));
      REQUIRE(dict.nodes().count("list"));
    } //
    SECTION("mesh") {
      PolyMesh pm;
      pm.readBoundary(asset_path);
      REQUIRE(pm.patches.size() == 5);
      {
        std::string patch_name = "frontAndBack";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 2);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 0);
      } //
      {
        std::string patch_name = "terrain";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 1);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 2);
      } //
      {
        std::string patch_name = "top";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 1);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 3);
      } //
      {
        std::string patch_name = "inlet";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 1);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 4);
      } //
      {
        std::string patch_name = "outlet";
        REQUIRE(pm.containsPatch(patch_name));
        auto patch = pm.patch(patch_name);
        REQUIRE(patch.size == 1);
        REQUIRE(patch.name == patch_name);
        REQUIRE(patch.start == 5);
      } //
    } //
  } //
}

TEST_CASE("owner dict") {
  hermes::Path assets_path(ASSETS_PATH);
  SECTION("owner1") {
    hermes::Log::addOptions(hermes::logging_options::abbreviate);
    auto asset_path = assets_path + "owner1";
    SECTION("dict") {
      OpenFoamDict dict(asset_path);
      HERMES_LOG_VARIABLE(dict.nodes().size());
      for (const auto &node : dict.nodes()) {
        HERMES_LOG_VARIABLE(node.first);
      }
      REQUIRE(dict.nodes().count("FoamFile"));
      REQUIRE(dict.nodes().count("list"));
    } //
    SECTION("mesh") {
      PolyMesh pm;
      pm.readOwners(asset_path);
      REQUIRE(pm.owner.size() == 6);
      for (auto o : pm.owner)
        REQUIRE(o == 0);
    } //
  } //
  SECTION("owner2") {
    hermes::Log::addOptions(hermes::logging_options::abbreviate);
    auto asset_path = assets_path + "owner2";
    SECTION("dict") {
      OpenFoamDict dict(asset_path);
      HERMES_LOG_VARIABLE(dict.nodes().size());
      for (const auto &node : dict.nodes()) {
        HERMES_LOG_VARIABLE(node.first);
      }
      REQUIRE(dict.nodes().count("FoamFile"));
      REQUIRE(dict.nodes().count("list"));
    } //
    SECTION("mesh") {
      PolyMesh pm;
      pm.readOwners(asset_path);
      REQUIRE(pm.owner.size() == 6);
      for (size_t i = 1; i <= 6; ++i)
        REQUIRE(pm.owner[i - 1] == i);
    } //
  } //
}

TEST_CASE("scalar field") {
  hermes::Path assets_path(ASSETS_PATH);
  SECTION("fa_h") {
    hermes::Log::addOptions(hermes::logging_options::abbreviate);
    auto asset_path = assets_path + "fa_h";
    SECTION("dict") {
      OpenFoamDict dict(asset_path);
      dict.print();
      REQUIRE(dict.nodes().count("dimensions"));
      REQUIRE(dict.nodes().count("internalField"));
      REQUIRE(dict.nodes().count("boundaryField"));
      REQUIRE(dict["boundaryField"].fields.count("wall"));
      REQUIRE(dict["boundaryField"].fields.count("terrain"));
      REQUIRE(dict["boundaryField"].fields.count("top"));
      std::vector<double> data = OpenFoamDict::parseValuesFrom<double>(
          dict["boundaryField"]["terrain"]["value"].value);
      REQUIRE(data.size() == 10);
      for (int i = 1; i <= 10; ++i)
        REQUIRE_THAT(data[i - 1], Catch::Matchers::WithinAbs(i * 1e-5, 1e-8));
    } //
  } //
}
