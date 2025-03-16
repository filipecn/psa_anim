/// Copyright (c) 2023, FilipeCN.
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
#include <hermes/common/file_system.h>
#include <hermes/common/index.h>
#include <hermes/logging/logging.h>
#include <hermes/numeric/numeric.h>
#include <psa_anim/interpolators.h>
#include <psa_anim/partition_grid.h>

#include <interpolation.h>
#include <mutex>
#include <rasterizer.h>
#include <utils.h>

using namespace psa_anim;

TEST_CASE("simple", "[interpolation]") {

  Poly func;

  hermes::Str p_, f_, g_, P_, F_;

  std::vector<double> values;
  std::vector<hermes::vec3> grad_values;
  std::vector<hermes::point3> centers;
  std::vector<size_t> ids;
  {
    size_t neval = 1;
    double dx = 1.0; // 1. / (neval * 10);
    for (auto ij : hermes::range3(hermes::size3(neval + 1))) {
      hermes::point3 p(ij.j * dx, ij.i * dx, ij.k * dx);
      auto fp = func.f(p);
      auto gp = func.g(p);
      centers.emplace_back(p);
      values.emplace_back(fp);
      grad_values.emplace_back(gp.x, gp.y, gp.z);
      ids.emplace_back(ids.size());
      p_.appendLine(p.x, ", ", p.y, ", ", p.z);
      f_.appendLine(fp);
      g_.appendLine(gp.x, ", ", gp.y, ", ", gp.z);
    }
  }

  std::vector<hermes::point3> samples;
  std::vector<double> analytic;
  {
    size_t neval = 8;
    double dx = 1. / neval;
    for (auto ij : hermes::range3(hermes::size3(neval + 1))) {
      // hermes::point3 p(ij.i * dx, ij.j * dx, ij.k * dx);
      hermes::point3 p(0.5, 0.5, 0.5);
      samples.emplace_back(p);
      analytic.emplace_back(func.f(p));
      P_.appendLine(p.x, ", ", p.y, ", ", p.z);
      F_.appendLine(func.f(p));
      break;
    }
  }

  hermes::FileSystem::writeFile("centers", p_.str());
  hermes::FileSystem::writeFile("centers_f", f_.str());
  hermes::FileSystem::writeFile("centers_g", g_.str());
  hermes::FileSystem::writeFile("points", P_.str());
  hermes::FileSystem::writeFile("points_f", F_.str());

  // SECTION("MultiQuadricRBF") {
  //   psa_anim::InterpolationSystem<psa_anim::MultiQuadricRBF> system;

  //  system.init(centers, ids, 6);
  //  system.build(values, grad_values, ids);
  //  auto cvalues = system.eval(centers, ids, centers, 6);
  //  auto svalues = system.eval(centers, ids, samples, 6);

  //  for(size_t i = 0; i < centers.size(); ++i) {
  //    REQUIRE_THAT(values[i] - cvalues[i], Catch::Matchers::WithinAbs(0,
  //    1e-5)); REQUIRE_THAT(svalues[i] - analytic[i],
  //    Catch::Matchers::WithinAbs(0, 1e-5));
  //  }
  //}
  SECTION("CubicRBF") {
    psa_anim::InterpolationSystem<psa_anim::CubicRBF> system;

    system.init(centers, ids);
    system.build(values, grad_values, ids);
    auto cvalues = system.eval(centers, ids, centers);
    auto svalues = system.eval(centers, ids, samples);

    float error = 0;
    for (size_t i = 0; i < samples.size(); ++i) {
      std::cerr << svalues[i] << " != " << analytic[i] << std::endl;
      error = std::max(error, (float)std::fabs(svalues[i] - analytic[i]));
    }
    REQUIRE_THAT(error, Catch::Matchers::WithinAbs(0, 1e-5));
  }
}

TEST_CASE("fasshauer", "[interpolation]") {
  return;

  auto f = [](const hermes::point3 &p) -> double {
    return (std::tanh(9 * (p.y - p.x)) + 1) / (std::tanh(9) + 1);
  };
  auto fx = [](const hermes::point3 &p) -> double {
    auto tn = std::tanh(9 * (p.y - p.x));
    return 9 * (tn * tn - 1) / (std::tan(9) + 1);
  };
  auto fy = [](const hermes::point3 &p) -> double {
    auto tn = std::tanh(9 * (p.y - p.x));
    return 9 * (1 - tn * tn) / (std::tan(9) + 1);
  };

  std::vector<double> values;
  std::vector<hermes::vec3> grad_values;
  std::vector<hermes::point3> centers;
  std::vector<size_t> ids;
  {
    size_t neval = 4;
    double dx = 1. / neval;
    for (auto ij : hermes::range2(hermes::size2(neval + 1))) {
      hermes::point3 p(ij.i * dx, ij.j * dx, 0);
      centers.emplace_back(p);
      values.emplace_back(f(p));
      grad_values.emplace_back(fx(p), fy(p), 0);
      ids.emplace_back(ids.size());
    }
  }

  std::vector<hermes::point3> samples;
  std::vector<double> analytic;
  {
    size_t neval = 8;
    double dx = 1. / neval;
    for (auto ij : hermes::range2(hermes::size2(neval + 1))) {
      hermes::point3 p(ij.i * dx, ij.j * dx, 0);
      samples.emplace_back(p);
      analytic.emplace_back(f(p));
    }
  }

  SECTION("MultiQuadricRBF") {
    psa_anim::InterpolationSystem<psa_anim::MultiQuadricRBF> system;

    system.init(centers, ids, 6);
    system.build(values, grad_values, ids);
    auto cvalues = system.eval(centers, ids, centers, 6);
    auto svalues = system.eval(centers, ids, samples, 6);

    for (size_t i = 0; i < centers.size(); ++i) {
      REQUIRE_THAT(values[i] - cvalues[i], Catch::Matchers::WithinAbs(0, 1e-5));
      REQUIRE_THAT(svalues[i] - analytic[i],
                   Catch::Matchers::WithinAbs(0, 1e-5));
    }
  }
  SECTION("CubicRBF") {
    psa_anim::InterpolationSystem<psa_anim::CubicRBF> system;

    system.init(centers, ids);
    system.build(values, grad_values, ids);
    auto cvalues = system.eval(centers, ids, centers);
    auto svalues = system.eval(centers, ids, samples);

    for (size_t i = 0; i < centers.size(); ++i) {
      REQUIRE_THAT(values[i] - cvalues[i], Catch::Matchers::WithinAbs(0, 1e-5));
      REQUIRE_THAT(svalues[i] - analytic[i],
                   Catch::Matchers::WithinAbs(0, 1e-5));
    }
  }
}

TEST_CASE("Cell Interpolator", "[interpolation]") {
  // setup mesh
  hermes::size3 res(20, 20, 20);
  hermes::vec3 mesh_dx(2);
  auto r_mesh = psa_anim::create_openfoam_grid_mesh(res, mesh_dx, 90 /*angle*/);
  r_mesh->computeCellBounds();
  psa_anim::PartitionGrid partition(*r_mesh, 10);

  CellInterpolator interpolator;

  std::vector<std::pair<size_t, size_t>> cells = {{0, 2}, {1, 5}, {2, 7}};
  std::vector<double> values = {0, 1, 2};
  interpolator.update(partition, {}, values, {});
  auto ivalues = interpolator.eval(partition, {}, {}, cells);
  HERMES_LOG_ARRAY(ivalues);
}

TEST_CASE("parallel", "[interpolation]") {
  std::mutex m;
  for_indices(100, 10, [&](size_t start, size_t end) {
    std::unique_lock<std::mutex> lock(m);
    REQUIRE(end - start + 1 == 10);
  });
}

TEST_CASE("partition", "[interpolation]") {
  // setup mesh
  hermes::size3 res(20, 20, 20);
  hermes::vec3 mesh_dx(2);
  auto r_mesh = psa_anim::create_openfoam_grid_mesh(res, mesh_dx, 90 /*angle*/);
  r_mesh->computeCellBounds();
  SECTION("transform") {
    f32 cell_size = 2;
    hermes::vec3 offset(2, 0, 0);
    f32 rcp = 1 / cell_size;
    auto toGrid = hermes::Transform::scale(rcp, rcp, rcp) *
                  hermes::Transform::translate(-offset);
    auto toWorld = hermes::Transform::translate(offset) *
                   hermes::Transform::scale(cell_size, cell_size, cell_size);
    HERMES_LOG_VARIABLE(toWorld(toGrid(hermes::point3(0, 0, 0))));
    HERMES_LOG_VARIABLE(toWorld(toGrid(hermes::point3(-1, 0, 0))));
    HERMES_LOG_VARIABLE(toWorld(toGrid(hermes::point3(1, 0, 0))));
    HERMES_LOG_VARIABLE(toWorld(toGrid(hermes::point3(0, -1, 0))));
    HERMES_LOG_VARIABLE(toWorld(toGrid(hermes::point3(0, 1, 0))));
    HERMES_LOG_VARIABLE(toWorld(toGrid(hermes::point3(0, 0, -1))));
    HERMES_LOG_VARIABLE(toWorld(toGrid(hermes::point3(0, 0, 1))));
    auto gp = toGrid(hermes::point3(10, 10, 10));
    hermes::index3 ijk(gp.x < 0 ? -hermes::Numbers::ceil2Int(-gp.x) : (int)gp.x,
                       gp.y < 0 ? -hermes::Numbers::ceil2Int(-gp.y) : (int)gp.y,
                       gp.z < 0 ? -hermes::Numbers::ceil2Int(-gp.z)
                                : (int)gp.z);
    HERMES_LOG_VARIABLES(gp, ijk);
    exit(0);
  }
  SECTION("offset 0") {
    psa_anim::PartitionGrid partition(*r_mesh, 10);
    HERMES_LOG_VARIABLE(r_mesh->cell_centers[0]);
    auto regions = partition.init({0});
    HERMES_LOG_VARIABLE(partition.regionBox(0).lower);
  }
  SECTION("offset x") {
    psa_anim::PartitionGrid partition(*r_mesh, 2, {2.0, 0.0, 0.0});
    HERMES_LOG_VARIABLE(r_mesh->cell_centers[0]);
    auto regions = partition.init({0});
  }
  SECTION("sanity") {
    // setup partitions
    psa_anim::PartitionGrid partition(*r_mesh, 5);
    // get all cell ids
    std::vector<size_t> cell_ids;
    for (size_t cell_id = 0; cell_id < r_mesh->cells.size(); ++cell_id)
      cell_ids.emplace_back(cell_id);
    // setup
    auto regions = partition.init(cell_ids);
    REQUIRE(regions.size() == 512);
    for (auto region : regions) {
      auto cells = partition.regionCells(region);
    }
  }
}

TEST_CASE("interpolators", "[interpolation]") {
  // setup mesh
  hermes::size3 res(20, 20, 20);
  hermes::vec3 mesh_dx(2);
  auto r_mesh = psa_anim::create_openfoam_grid_mesh(res, mesh_dx, 90 /*angle*/);
  r_mesh->computeCellBounds();
  // setup partitions
  psa_anim::PartitionGrid partition(*r_mesh, 5);
  // get all cell ids
  std::vector<size_t> cell_ids;
  for (size_t cell_id = 0; cell_id < r_mesh->cells.size(); ++cell_id)
    cell_ids.emplace_back(cell_id);
  auto regions = partition.init(cell_ids);
  // fill
  Gaussian gaussian;
  auto values = gaussian.F(r_mesh->cell_centers);
  auto grad_values = gaussian.G(r_mesh->cell_centers);
  SECTION("RBF") {
    psa_anim::RBFInterpolator interpolator;
    interpolator.prepare(partition, regions);
    interpolator.update(partition, regions, values, grad_values);
    double rmse = 0.;
    for (auto region : regions) {
      auto p = partition.regionBox(region).center();
      auto vs = interpolator.eval(partition, region, {p});
      auto a = gaussian.f(p);
      rmse += (a - vs.front()) * (a - vs.front());
    }
    rmse = std::sqrt(rmse / regions.size());
    HERMES_LOG_VARIABLE(rmse);
  }
}

TEST_CASE("convergence", "[interpolation]") {
  return;
  std::vector<size_t> resolutions = {16};
  for (auto r : resolutions) {
    hermes::size3 res(r, r, r);
    hermes::vec3 mesh_dx(1. / r);
    auto r_mesh = psa_anim::create_openfoam_grid_mesh(res, mesh_dx, 90);
    r_mesh->computeCellBounds();
    // setup partitions
    psa_anim::PartitionGrid partition(*r_mesh, 5);
    // get all cell ids
    std::vector<size_t> cell_ids;
    for (size_t cell_id = 0; cell_id < r_mesh->cells.size(); ++cell_id)
      cell_ids.emplace_back(cell_id);
    auto regions = partition.init(cell_ids);
  }
}

TEST_CASE("raster", "[interpolation]") {
  hermes::size3 res(20, 20, 20);
  hermes::vec3 mesh_dx(2);
  auto r_mesh = psa_anim::create_openfoam_grid_mesh(res, mesh_dx, 90 /*angle*/);
  REQUIRE(r_mesh);
  REQUIRE(r_mesh->cells.size() == res.total());
  r_mesh->computeCellBounds();
  Gaussian gaussian;
  auto values = gaussian.F(r_mesh->cell_centers);
  auto grad_values = gaussian.G(r_mesh->cell_centers);
  psa_anim::Rasterizer<psa_anim::RBFInterpolator<psa_anim::CubicRBF>> raster(
      *r_mesh, 1, 5);
  std::vector<size_t> cell_ids;
  for (size_t i = 0; i < r_mesh->cells.size(); i += 100)
    cell_ids.emplace_back(i);
  raster.setup(cell_ids);

  std::mutex m;
  double sum = 0;
  raster.raster(cell_ids, values, grad_values,
                [&](const psa_anim::VoxelCallbackData &data) {
                  for (const auto &item : data.regions) {
                    double rmse = 0.0;
                    for (size_t i = 0; i < item.second.points.size(); ++i) {
                      auto error = item.second.values[i] -
                                   gaussian.f(item.second.points[i]);
                      rmse += error * error;
                    }
                    rmse = std::sqrt(rmse / item.second.points.size());
                    std::lock_guard<std::mutex> lock(m);
                    HERMES_LOG_ERROR("{} -> {} {}", item.first, rmse,
                                     item.second.points.size());
                  }
                });
  // double rmse = std::sqrt(sum / n);
  // HERMES_LOG_VARIABLE(rmse);
}
