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

#include <psa_anim/block_mesh.h>
#include <psa_anim/mesh_utils.h>
#include <psa_anim/openfoam_poly_mesh.h>

hermes::Path assets_path(ASSETS_PATH);

TEST_CASE("mesh2grid2") {
  auto terrain = psa_anim::loadSTL(
      "/mnt/windows/Projects/PHD/results/wolfsgrube2/dsl/constant/surface.stl",
      "terrain");
  REQUIRE(terrain);
  hermes::point2 a(1755, 40);
  hermes::point2 b(1350, 780);
  hermes::vec2 axis = hermes::normalize(b - a);
  hermes::point2 origin = a; // + axis * 100;
                             //  hermes::vec2 axis(-1162, 2456);
  auto grid =
      psa_anim::extractGrid(*terrain, axis, a, {(b - origin).length(), 600}, 5);
  psa_anim::info(grid);
  psa_anim::saveOBJ(grid, "/home/filipecn/Desktop/grid.obj");
}

TEST_CASE("mesh2grid") {
  auto terrain =
      psa_anim::loadSTL("/mnt/windows/Projects/foam_experiments/avalanche/"
                        "tutorials/wolfsgrube/constant/surface.stl",
                        "terrain");
  REQUIRE(terrain);
  hermes::point2 a(1650, 150);
  hermes::point2 b(600, 2040);
  hermes::vec2 axis = hermes::normalize(b - a);
  hermes::point2 origin = a; // + axis * 100;
                             //  hermes::vec2 axis(-1162, 2456);
  auto grid = psa_anim::extractGrid(*terrain, axis, a,
                                    {(b - origin).length() * 0.95f, 600}, 5);
  psa_anim::info(grid);
  psa_anim::saveOBJ(grid, "/home/filipecn/Desktop/grid.obj");
}

TEST_CASE("block_mesh") {
  SECTION("sanity") {
    // boundaries
    std::map<std::string, std::vector<hermes::vec3>> boundaries = {
        {"top", {{0, 0, 1}}},
        {"bottom", {{0, 0, -1}}},
        {"outlet", {{1, 0, 0}}},
        {"inlet", {{-1, 0, 0}}},
        {"walls", {{0, 1, 0}, {0, -1, 0}}}};

    auto bm = psa_anim::blockMesh(
        {
            {0, 0, 0}, // vertex 0
            {1, 0, 0}, // vertex 1
            {2, 0, 0}, // vertex 2
            {0, 1, 0}, // vertex 3
            {1, 1, 0}, // vertex 4
            {2, 1, 0}  // vertex 5
        },
        {
            {0, 3, 4, 1}, // face 0
            {1, 4, 5, 2}  // face 1
        },
        boundaries);

    bm.save("/mnt/windows/Projects/PHD/results/wolfsgrube2/psl/system/"
            "blockMeshDict");
  } //
  SECTION("grid") {
    // boundaries
    std::map<std::string, std::vector<hermes::vec3>> boundaries = {
        {"top", {{0, 0, 1}}},
        {"bottom", {{0, 0, -1}}},
        {"outlet", {{1, 0, 0}}},
        {"inlet", {{-1, 0, 0}}},
        {"walls", {{0, 1, 0}, {0, -1, 0}}}};

    std::vector<hermes::point3> points;
    std::vector<std::vector<size_t>> faces;

    auto grid_size = 5;

    for (auto ij : hermes::range2(hermes::size2(grid_size, grid_size))) {
      points.emplace_back(ij.i, ij.j, 0);
    }
    auto indexOf = [](const hermes::Index2<i32> &ij, size_t w) {
      return ij.j * w + ij.i;
    };
    for (auto ij :
         hermes::range2(hermes::size2(grid_size - 1, grid_size - 1))) {
      faces.push_back({
          indexOf(ij.plus(0, 0), grid_size),
          indexOf(ij.plus(0, 1), grid_size),
          indexOf(ij.plus(1, 1), grid_size),
          indexOf(ij.plus(1, 0), grid_size),
      });
    }

    auto bm = psa_anim::blockMesh(points, faces, boundaries, 10);

    bm.save("/mnt/windows/Projects/PHD/results/wolfsgrube2/psl/system/"
            "blockMeshDict");
  } //
  SECTION("procedural terrain") {
    // boundaries
    std::map<std::string, std::vector<hermes::vec3>> boundaries = {
        {"top", {{0, 0, 1}}},
        {"bottom", {{0, 0, -1}}},
        {"outlet", {{1, 0, 0}}},
        {"inlet", {{-1, 0, 0}}},
        {"walls", {{0, 1, 0}, {0, -1, 0}}}};

    std::vector<hermes::point3> points;
    std::vector<std::vector<size_t>> faces;

    auto grid_size = 5;

    for (auto ij : hermes::range2(hermes::size2(grid_size, grid_size))) {
      points.emplace_back(ij.i, ij.j, 0);
    }
    auto indexOf = [](const hermes::Index2<i32> &ij, size_t w) {
      return ij.j * w + ij.i;
    };
    for (auto ij :
         hermes::range2(hermes::size2(grid_size - 1, grid_size - 1))) {
      faces.push_back({
          indexOf(ij.plus(0, 0), grid_size),
          indexOf(ij.plus(0, 1), grid_size),
          indexOf(ij.plus(1, 1), grid_size),
          indexOf(ij.plus(1, 0), grid_size),
      });
    }

    auto bm = psa_anim::blockMesh(points, faces, boundaries, 10);

    bm.save("/mnt/windows/Projects/PHD/results/wolfsgrube2/psl/system/"
            "blockMeshDict");
  } //
  SECTION("from obj") {
    //    auto mesh = psa_anim::loadOBJ("/home/filipecn/Desktop/quadm.obj");
    auto mesh = psa_anim::loadOBJ("/home/filipecn/Desktop/test.obj");
    REQUIRE(mesh);
    psa_anim::info(*mesh);
    std::map<std::string, std::vector<hermes::vec3>> boundaries = {
        {"top", {{0, 0, 1}}},
        {"bottom", {{0, 0, -1}}},
        {"outlet", {{1, 0, 0}}},
        {"inlet", {{-1, 0, 0}}},
        {"walls", {{0, 1, 0}, {0, -1, 0}}}};
    auto bm = psa_anim::blockMesh(mesh->vertices, mesh->faces, boundaries, 10);
    bm.save("/mnt/windows/Projects/PHD/results/wolfsgrube2/psl/system/"
            "blockMeshDict");
  } //
  SECTION("from stl") {
    auto mesh = psa_anim::loadSTL("/home/filipecn/Desktop/test.stl");
    REQUIRE(mesh);
    //    auto quad = psa_anim::convertToQuad(*mesh);
  } //
}
