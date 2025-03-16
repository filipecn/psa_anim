/// Copyright (c) 2022, FilipeCN.
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
///\file foam2vdb.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-03-11
///
///\brief

#include <hermes/common/arg_parser.h>
#include <hermes/common/file_system.h>
#include <hermes/common/profiler.h>
#include <hermes/geometry/bbox.h>
#include <hermes/logging/logging.h>
#include <interpolators.h>
#include <mesh_utils.h>
#include <mutex>
#include <openvdb/openvdb.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <partition_grid.h>
#include <rasterizer.h>
#include <sim_path.h>
#include <unistd.h>
#include <utils.h>
#include <volume.h>

using namespace psa_anim;
using RBFVolume = Volume<Rasterizer<RBFInterpolator<CubicRBF>>>;
using HRBFVolume = Volume<Rasterizer<HRBFInterpolator<CubicRBF>>>;
using CellVolume = Volume<Rasterizer<CellInterpolator>>;
using ShepardVolume = Volume<Rasterizer<ShepardInterpolator>>;
using NearestVolume = Volume<Rasterizer<NearestInterpolator>>;

enum class InterpolationMethod { HRBF, CELL, SHEPARD, RBF, NEAREST };

std::string methodName(InterpolationMethod method) {
  if (method == InterpolationMethod::RBF)
    return "rbf";
  if (method == InterpolationMethod::HRBF)
    return "hrbf";
  if (method == InterpolationMethod::SHEPARD)
    return "shepard";
  if (method == InterpolationMethod::NEAREST)
    return "nearest";
  return "cell";
}

struct Parameters {
  // input/output config
  hermes::Path path;
  hermes::Path output_path;
  std::string grid_name;
  int single_frame;
  bool force_write;
  bool only_missing_frames;
  double dt;
  bool renumber_frames;
  bool renumber_frames_seq;
  bool smooth_data_field;
  int every_n;
  bool gen_isosurface;
  float is_2d;
  bool only_cells;
  // conversion
  double isovalue;
  double min_cut;
  double voxel_size;
  // debug
  bool debug;
  // interpolation
  InterpolationMethod interpolation_method{InterpolationMethod::HRBF};
  // ptu
  bool use_ptu;
  bool ptu_diag;
  float ptu_size;
  // test
  bool test;
  bool test_section;
  float test_angle;
  int test_res;
  float test_dx;
  float test_dy;
  float test_gaussian_c;
  float test_gaussian_a;
  bool test_animate_scale;
  bool test_animate_angle;
  std::string color_by;

  int read(int argc, const char **argv) {
    hermes::ArgParser parser("foam2vdb",
                             "converts a OpenFOAM volume field simulation data "
                             "into a OpenVDB volume animation frames.");
    parser.addArgument("-i", "OpenFOAM simulation path");
    parser.addArgument("-o", "OpenVDB output path");
    parser.addArgument("-n", "vdb grid name [default = density]");
    parser.addArgument("--dt",
                       "Re-number frames based on a given dt [default = 0.05]");
    parser.addArgument("--renumber-frames",
                       "Re-number frames based on the given dt");
    parser.addArgument("--renumber-frames-sequentially",
                       "Re-number frames sequentially");
    parser.addArgument("--color-by", "[density | velocity]");
    parser.addArgument("--smooth", "Smooth input data field");
    parser.addArgument("-f", "single frame number");
    parser.addArgument("-F", "force write");
    parser.addArgument("--only-missing-frames", "process only missing files");
    parser.addArgument("--every-n", "convert every n-th");
    parser.addArgument("--isosurface", "gen surface");
    parser.addArgument("-2d", "is 2d");
    parser.addArgument("--only-cells", "gen voxels only in density cells");
    // conversion
    parser.addArgument("-s", "OpenVDB voxel size [default = 0.2]");
    parser.addArgument("--min", "Filter values bellow this [default = 0.001]");
    parser.addArgument("--iso", "Layer interface isovalue [default = 0.3]");
    // debug
    parser.addArgument("--debug", "gen debug data");
    // interpolation
    parser.addArgument("--interpolation", "method: cell or rbf");
    // ptu
    parser.addArgument("--use-ptu", "use partition of unity");
    parser.addArgument("--ptu-size", "partition of unity region size");
    parser.addArgument("--ptu-diag", "also use diagonal shift in ptu");
    // test
    parser.addArgument("--test", "run test");
    parser.addArgument("--test-section", "run test section");
    parser.addArgument("--test-res", "run test res");
    parser.addArgument("--test-angle", "run test angle");
    parser.addArgument("--test-dx", "test cell size in x");
    parser.addArgument("--test-dy", "test cell size in z");
    parser.addArgument("--test-gaussian-c", "test gaussian c");
    parser.addArgument("--test-gaussian-a", "test gaussian a");
    parser.addArgument("--test-animate-angle", "test varying mesh angle");
    parser.addArgument("--test-animate-scale", "test varying mesh angle");
    parser.parse(argc, argv, true);
    // parameters
    path = parser.get<std::string>("-i");
    output_path = parser.get<std::string>("-o");
    grid_name = parser.get<std::string>("-n", "density");
    single_frame = parser.get<int>("-f", -1);
    force_write = parser.check("-F");
    only_missing_frames = parser.check("--only-missing-frames");
    dt = parser.get<float>("--dt", 0.05);
    renumber_frames = parser.check("--renumber-frames");
    renumber_frames_seq = parser.check("--renumber-frames-sequentially");
    smooth_data_field = parser.check("--smooth");
    every_n = parser.get<int>("--every-n", 1);
    gen_isosurface = parser.check("--isosurface");
    color_by = parser.get<std::string>("--color-by", "");
    is_2d = parser.check("-2d");
    only_cells = parser.check("--only-cells");
    // debug
    debug = parser.check("--debug");
    // conversion
    voxel_size = parser.get<float>("-s", 0.5);
    min_cut = parser.get<float>("--min", 0.001);
    isovalue = parser.get<float>("--iso", 0.001);
    // interpolation
    std::string interpolation =
        parser.get<std::string>("--interpolation", "hrbf");
    // ptu
    use_ptu = parser.check("--use-ptu");
    ptu_diag = parser.check("--ptu-diag");
    ptu_size = parser.get<float>("--ptu-size", 20);
    // test
    test = parser.check("--test");
    test_section = parser.check("--test-section");
    test_angle = parser.get<float>("--test-angle", 90);
    test_res = parser.get<int>("--test-res", 20);
    test_dx = parser.get<float>("--test-dx", 1.f);
    test_dy = parser.get<float>("--test-dy", 1.f);
    test_gaussian_c = parser.get<float>("--test-gaussian-c", 5.0);
    test_gaussian_a = parser.get<float>("--test-gaussian-a", 1.0);
    test_animate_angle = parser.check("--test-animate-angle");
    test_animate_scale = parser.check("--test-animate-scale");

    // check output directory
    if (!test && !output_path.isDirectory()) {
      HERMES_LOG_ERROR("invalid output path");
      return 0;
    }
    if (!test && !test_section && !only_missing_frames && single_frame < 0 &&
        !hermes::FileSystem::ls(output_path).empty()) {
      HERMES_LOG_ERROR("output directory not empty");
      return 0;
    }

    if (interpolation == "rbf")
      interpolation_method = InterpolationMethod::RBF;
    else if (interpolation == "hrbf")
      interpolation_method = InterpolationMethod::HRBF;
    else if (interpolation == "cell")
      interpolation_method = InterpolationMethod::CELL;
    else if (interpolation == "shepard")
      interpolation_method = InterpolationMethod::SHEPARD;
    else if (interpolation == "nearest")
      interpolation_method = InterpolationMethod::NEAREST;
    else {
      HERMES_LOG_ERROR("invalid interpolation method");
      return 0;
    }
    return 1;
  }
} args;

int volume_count = 5;

// debug
struct DebugData {
  void clear() {
    region_ids.clear();
    voxels.clear();
    voxel_values.clear();
  }
  std::set<size_t> region_ids;
  std::unordered_map<size_t, std::vector<hermes::point3>> voxels;
  std::unordered_map<size_t, std::vector<double>> voxel_values;
};

struct VDBVolumeInterface {
  virtual void update(size_t frame, const std::vector<double> &values,
                      const std::vector<hermes::vec3> &grad_values) = 0;

  virtual openvdb::FloatGrid::Ptr output(size_t frame,
                                         const hermes::Path &path) = 0;

  virtual void dumpDebug(size_t frame) = 0;

  std::vector<DebugData> debug_data;
};

template <typename VolumeType> struct PUVolume : public VDBVolumeInterface {
  PUVolume(const PolyMesh &mesh) : mesh(mesh) {
    a_palette = psa_anim::ColorPalettes::Batlow();
    min_a = 0.0;
    max_a = 1.0;
    HERMES_PROFILE_FUNCTION();
    f32 s_2 = args.ptu_size * 0.5;
    volumes.emplace_back(
        new VolumeType(mesh, args.voxel_size, args.ptu_size, {0.0, 0.0, 0.0}));
    if (args.use_ptu) {
      volumes.emplace_back(new VolumeType(mesh, args.voxel_size, args.ptu_size,
                                          {s_2, 0.0, 0.0}));
      volumes.emplace_back(new VolumeType(mesh, args.voxel_size, args.ptu_size,
                                          {0.0, s_2, 0.0}));
      volumes.emplace_back(new VolumeType(mesh, args.voxel_size, args.ptu_size,
                                          {0.0, 0.0, s_2}));
      if (args.ptu_diag)
        volumes.emplace_back(new VolumeType(mesh, args.voxel_size,
                                            args.ptu_size, {s_2, s_2, s_2}));
    }
    volume_count = volumes.size();
    RAM();
    HERMES_LOG("volume count {}", volume_count);
    if (args.debug)
      debug_data.resize(volume_count);
  }

  void update(size_t frame, const std::vector<double> &values,
              const std::vector<hermes::vec3> &grad_values) override {
    HERMES_PROFILE_FUNCTION();
    ThreadPool pool(volumes.size());
    for (size_t i = 0; i < volumes.size(); ++i) {
      pool.enqueue([i, this, &frame, &values, &grad_values]() {
        volumes[i]->rasterizer().set2d(args.is_2d);
        volumes[i]->rasterizer().setOnlyCells(args.only_cells);
        std::vector<size_t> pu_regions;
        if (args.debug) {
          debug_data[i].clear();
          openvdb::FloatGrid::Ptr vdb = openvdb::FloatGrid::create();
          vdb->setTransform(volumes[i]->transform());
          auto acc = vdb->getAccessor();

          std::mutex voxel_mutex;
          pu_regions = volumes[i]->fill(
              frame, args.min_cut, values, grad_values,
              [&](size_t region_id, const hermes::point3 &p,
                  double &v) -> bool {
                std::lock_guard<std::mutex> lock(voxel_mutex);
                debug_data[i].region_ids.insert(region_id);
                debug_data[i].voxels[region_id].emplace_back(p);
                debug_data[i].voxel_values[region_id].emplace_back(v);
                acc.setValue(openvdb::Coord(p.x, p.y, p.z), v);
                return false;
              });
          openvdb::io::File(hermes::Str::concat("volume", i, ".vdb"))
              .write({vdb});

        } else {
          pu_regions =
              volumes[i]->fill(frame, args.min_cut, values, grad_values);
        }
        HERMES_LOG_VARIABLE(pu_regions.size());
        size_t mean_cell_count = 0;
        size_t min_cell_count = 1 << 20;
        size_t max_cell_count = 0;
        for (auto pu_region : pu_regions) {
          size_t cell_count = volumes[i]
                                  ->rasterizer()
                                  .interpolationGrid()
                                  .regionCells(pu_region)
                                  .size();
          min_cell_count = std::min(min_cell_count, cell_count);
          max_cell_count = std::max(max_cell_count, cell_count);
          mean_cell_count += cell_count;
        }
        if (!pu_regions.empty())
          mean_cell_count /= pu_regions.size();
        HERMES_LOG_VARIABLE(mean_cell_count);
        HERMES_LOG_VARIABLES(min_cell_count, max_cell_count);
        auto face_mesh =
            volumes[i]->rasterizer().interpolationGrid().regionCellsMesh(0);
        psa_anim::saveOBJ(face_mesh, hermes::Str::concat("cells", i, ".obj"),
                          true);
        psa_anim::bbox2obj(
            {volumes[i]->rasterizer().interpolationGrid().regionBox(0)},
            hermes::Str::concat("grid", i, ".obj"));
      });
    }
  }

  struct Local {
    static inline void write(openvdb::CombineArgs<float> &args) {
      // Transfer the B value and its active state.
      args.setResult(args.b() / volume_count);
      args.setResultIsActive(args.bIsActive());
    }
    static inline void blend(openvdb::CombineArgs<float> &args) {
      args.setResult(args.a() + args.b() / volume_count);
      args.setResultIsActive(args.bIsActive() || args.aIsActive());
    }
  };

  static void testBlend() {
    {
      volume_count = 2;
      // A and B
      openvdb::FloatGrid::Ptr A = openvdb::FloatGrid::create();
      openvdb::FloatGrid::Ptr B = openvdb::FloatGrid::create();
      openvdb::FloatGrid::Accessor a = A->getAccessor();
      openvdb::FloatGrid::Accessor b = B->getAccessor();
      a.setValue({0, 0, 0}, 4);
      a.setValue({1, 0, 0}, 4);
      a.setValue({2, 0, 0}, 4);
      b.setValue({1, 0, 0}, 8);
      b.setValue({2, 0, 0}, 8);
      b.setValue({3, 0, 0}, 8);
      openvdb::FloatGrid::Ptr vdb = openvdb::FloatGrid::create();
      vdb->tree().combineExtended(A->tree(), Local::write);
      // Iterate over all active values but don't allow them to be changed.
      for (openvdb::FloatGrid::ValueOnCIter iter = vdb->cbeginValueOn();
           iter.test(); ++iter) {
        if (iter.isVoxelValue()) {
          std::cout << iter.getCoord() << " " << *iter << std::endl;
        } else {
          openvdb::CoordBBox bbox;
          iter.getBoundingBox(bbox);
          std::cout << bbox << std::endl;
        }
      }

      vdb->tree().combineExtended(B->tree(), Local::blend);

      // Iterate over all active values but don't allow them to be changed.
      for (openvdb::FloatGrid::ValueOnCIter iter = vdb->cbeginValueOn();
           iter.test(); ++iter) {
        if (iter.isVoxelValue()) {
          std::cout << iter.getCoord() << " " << *iter << std::endl;
        } else {
          openvdb::CoordBBox bbox;
          iter.getBoundingBox(bbox);
          std::cout << bbox << std::endl;
        }
      }
    }
  }

  openvdb::FloatGrid::Ptr output(size_t frame,
                                 const hermes::Path &path) override {
    HERMES_LOG("building final vdb grid");
    openvdb::FloatGrid::Ptr vdb = openvdb::FloatGrid::create();
    vdb->setTransform(volumes[0]->transform());
    vdb->setName("density");
    RAM();
    for (size_t i = 0; i < volumes.size(); ++i) {
      auto v_vdb = volumes[i]->output(frame);
      if (i == 0)
        vdb->tree().combineExtended(v_vdb->tree(), Local::write);
      else
        vdb->tree().combineExtended(v_vdb->tree(), Local::blend);
    }
    HERMES_LOG("writing {}", path.fullName());
    RAM();
    // color
    if (args.color_by == "density") {
      openvdb::Vec3fGrid::Ptr color_vdb = openvdb::Vec3fGrid::create();
      color_vdb->setTransform(volumes[0]->transform());
      color_vdb->setName("color");
      openvdb::Vec3fGrid::Accessor a = color_vdb->getAccessor();
      double max_value = 0;
      double min_value = 1 << 20;
      size_t voxel_count = 0;
      for (openvdb::FloatGrid::ValueOnCIter iter = vdb->cbeginValueOn();
           iter.test(); ++iter) {
        if (iter.isVoxelValue()) {
          max_value = std::max(max_value, (double)*iter);
          min_value = std::min(min_value, (double)*iter);
          auto color = computeAColor(*iter);
          a.setValue(iter.getCoord(), {color.r, color.g, color.b});
          voxel_count++;
        }
      }
      HERMES_LOG_WARNING("DENSITY RANGE [{} {}] {} voxels", min_value,
                         max_value, voxel_count);
      openvdb::io::File(path.fullName()).write({vdb, color_vdb});
    } else {
      openvdb::io::File(path.fullName()).write({vdb});
    }
    return vdb;
    if (args.gen_isosurface && vdb->activeVoxelCount()) {
      HERMES_LOG("building isosurface");

      std::vector<openvdb::Vec3s> points;
      std::vector<openvdb::Vec4I> quads;
      std::vector<openvdb::Vec3I> triangles;

      openvdb::v11_0::tools::volumeToMesh(*vdb, points, triangles, quads,
                                          args.isovalue);
      // store surface
      hermes::Str surf_file_content;
      surf_file_content.appendLine("#v_f3_f4: ", points.size(), " ",
                                   triangles.size(), " ", quads.size());
      for (const auto &point : points)
        surf_file_content.appendLine("v ", point.x(), " ", point.y(), " ",
                                     point.z());
      for (const auto &triangle : triangles)
        surf_file_content.appendLine("f ", triangle.x() + 1, " ",
                                     triangle.y() + 1, " ", triangle.z() + 1);
      for (const auto &quad : quads)
        surf_file_content.appendLine("f ", quad.x() + 1, " ", quad.y() + 1, " ",
                                     quad.z() + 1, " ", quad.w() + 1);
      auto filename =
          hermes::FileSystem::basename(path.fullName(), ".vdb") + ".obj";
      HERMES_LOG("output file: [{}]", filename);
      hermes::FileSystem::writeFile(filename, surf_file_content.str());
    }
    return vdb;
  }

  void dumpDebug(size_t frame) override {
    HERMES_PROFILE_FUNCTION();
    static char volume_labels[5] = {'0', 'x', 'y', 'z', 'd'};
    auto face_mesh = polyMesh2faceMesh(mesh);
    saveOBJ(face_mesh, "mesh.obj", true);
    // points2ply(mesh.cell_centers, "cells.ply", values, false);
    ThreadPool pool(volumes.size());
    for (size_t i = 0; i < debug_data.size(); ++i) {
      pool.enqueue([i, this]() {
        auto prefix = hermes::Str::concat(volume_labels[i], "_");
        std::vector<hermes::point3> all_voxels;
        std::vector<double> all_values;
        for (const auto &item : debug_data[i].voxels)
          all_voxels.insert(all_voxels.end(), item.second.begin(),
                            item.second.end());
        for (const auto &item : debug_data[i].voxel_values)
          all_values.insert(all_values.end(), item.second.begin(),
                            item.second.end());
        std::vector<hermes::bbox3> regions;
        for (auto region_id : debug_data[i].region_ids)
          regions.emplace_back(volumes[i]->interpolationBox(region_id));
        std::unordered_map<size_t, std::vector<size_t>> stencils;
        for (auto id : debug_data[i].region_ids)
          stencils[id] =
              volumes[i]->rasterizer().interpolationGrid().regionCells(id);
        bbox2obj(regions, prefix + "regions.obj");
        // bbox2ply(volumes[i].regions, "regions.ply", rmses);
        points2ply(all_voxels, prefix + "voxels.ply", all_values, false);
        stencils2ply(mesh.cell_centers, stencils, prefix + "stencils.ply");
      });
    }
  }

  std::vector<std::unique_ptr<VolumeType>> volumes;
  const PolyMesh &mesh;
  psa_anim::ColorPalette a_palette;
  float min_a{hermes::Numbers::greatest_f32()};
  float max_a{hermes::Numbers::lowest_f32()};
  psa_anim::Color computeAColor(float value) {
    return a_palette((value - min_a) / (max_a - min_a));
  }
};

struct VecPUVolume {
  VecPUVolume(const PolyMesh &mesh) {
    for (int i = 0; i < 3; ++i)
      volumes[i].reset(new PUVolume<RBFVolume>(mesh));
  }

  void update(size_t frame, const std::vector<hermes::vec3> &u_values) {
    std::vector<double> values[3];
    for (const auto &u : u_values)
      for (int i = 0; i < 3; ++i)
        values[i].emplace_back(u[i]);
    for (int i = 0; i < 3; ++i) {
      volumes[i]->update(frame, values[i], {});
      auto grid = volumes[i]->output(frame, "xyz");
    }
  }
  // xyz
  std::unique_ptr<PUVolume<RBFVolume>> volumes[3];
};

using RegionError = std::pair<size_t, double>;

std::vector<RegionError> RMSE(const Gaussian &gaussian, DebugData &data) {
  std::vector<RegionError> errors;
  for (const auto &region_id : data.region_ids) {
    const auto &voxels = data.voxels[region_id];
    const auto &values = data.voxel_values[region_id];
    double sum = 0;
    HERMES_ASSERT(voxels.size() == values.size());
    HERMES_ASSERT(!voxels.empty());
    for (size_t i = 0; i < voxels.size(); ++i) {
      auto ans = gaussian.f(voxels[i]);
      auto num = values[i];
      sum += (ans - num) * (ans - num);
    }
    errors.push_back(RegionError(region_id, std::sqrt(sum / voxels.size())));
  }
  return errors;
}

int test() {
  {
    HERMES_PROFILE_SCOPE("test");
    HERMES_LOG("Building mesh");

    hermes::Str csv;
    csv.appendLine("method, angle, cell_size, rmse");

    Gaussian gaussian;
    gaussian.a = args.test_gaussian_a;
    gaussian.c = args.test_gaussian_c;
    // gaussian.center = {5, 5, 5};
    // gaussian.standard_deviation = args.test_gaussian_c;
    // gaussian.mean = {};

    f32 domain_size = args.test_res * args.test_dx;
    hermes::size3 domain_res(args.test_res, args.test_res, args.test_res);
    size_t fps = 24;
    float duration = 3;
    size_t n_frames = (args.test_animate_scale || args.test_animate_angle)
                          ? fps * duration
                          : 1;
    hermes::vec3 cell_size{args.test_dx, args.test_dx, args.test_dy};
    float angle = args.test_angle;

    for (size_t frame = 0; frame < n_frames; ++frame) {
      HERMES_LOG("Computing frame {}", frame);
      if (args.test_animate_angle)
        angle = 90 - (90 - 30) * (frame / (double)n_frames);
      if (args.test_animate_scale)
        cell_size.x = 1 + (10.0 * frame / (double)n_frames);
      HERMES_LOG_VARIABLES(cell_size, angle);

      // create grid
      auto r_mesh = create_openfoam_grid_mesh(domain_res, cell_size, angle);
      r_mesh->computeCellBounds();
      r_mesh->computeCellCenters();

      // setup vdb
      VDBVolumeInterface *volume = nullptr;
      switch (args.interpolation_method) {
      case InterpolationMethod::HRBF:
        volume = new PUVolume<HRBFVolume>(*r_mesh);
        break;
      case InterpolationMethod::NEAREST:
        volume = new PUVolume<NearestVolume>(*r_mesh);
        break;
      case InterpolationMethod::RBF:
        volume = new PUVolume<RBFVolume>(*r_mesh);
        break;
      case InterpolationMethod::CELL:
        volume = new PUVolume<CellVolume>(*r_mesh);
        break;
      case InterpolationMethod::SHEPARD:
        volume = new PUVolume<ShepardVolume>(*r_mesh);
        break;
      default:
        return -1;
      }

      auto values = gaussian.F(r_mesh->cell_centers);
      auto grad_values = gaussian.G(r_mesh->cell_centers);
      // fill volume
      volume->update(frame, values, grad_values);
      // debug data
      volume->dumpDebug(frame);
      // compute error
      /*
      if (reinterpret_cast<PUVolume<HRBFVolume> *>(volume)) {
        auto *rbf_volume = reinterpret_cast<PUVolume<HRBFVolume> *>(volume);
        if (!rbf_volume->debug_data.empty()) {
          double mean_rmse = 0;
          for (size_t volume_id = 0; volume_id < volume->debug_data.size();
               ++volume_id) {
            auto &debug_data = volume->debug_data[volume_id];
            auto errors = RMSE(gaussian, debug_data);
            std::sort(errors.begin(), errors.end(),
                      [](const RegionError &a, const RegionError &b) {
                        return a.second > b.second;
                      });
            std::vector<hermes::bbox3> region_boxes;
            for (auto error : errors) {
              region_boxes.push_back(rbf_volume->volumes[volume_id]
                                         ->rasterizer()
                                         .interpolationGrid()
                                         .regionBox(error.first));
            }
            psa_anim::bbox2obj(
                region_boxes,
                hermes::Str::concat("volume", volume_id, "_regions.obj"));
          }
        }
      }*/
      // compute output
      auto fn = hermes::Str::concat(
          // method
          methodName(args.interpolation_method), "_",
          // angle
          "angle_", std::setfill('0'), std::setw(6), (int)(angle * 1000), "_",
          // dx
          "dx_", std::setfill('0'), std::setw(6), (int)(cell_size.x * 1000),
          // dy
          "dy_", std::setfill('0'), std::setw(6), (int)(cell_size.y * 1000),
          "_",
          // frame
          std::setfill('0'), std::setw(6), frame);
      auto outputFile = hermes::Str::concat("test_", fn, ".vdb");
      if (!args.output_path.fullName().empty())
        outputFile =
            (args.output_path / hermes::Str::concat("test_", fn, ".vdb"))
                .fullName();
      auto vdb = volume->output(frame, outputFile);
      // compute final error
      double rmse = 0;
      size_t count = 0;
      double acc = 0;
      for (openvdb::FloatGrid::ValueOnIter iter = vdb->beginValueOn(); iter;
           ++iter) {
        if (iter.isVoxelValue()) {
          double value = *iter;
          acc += value;
          openvdb::Vec3f wp = vdb->transform().indexToWorld(iter.getCoord());
          // hermes::point3 p(iter.getCoord().x(), iter.getCoord().y(),
          //                  iter.getCoord().z());
          hermes::point3 p(wp.x(), wp.y(), wp.z());
          auto ans = gaussian.f(p);
          rmse += (value - ans) * (value - ans);
          count++;
        }
      }
      HERMES_ASSERT(count);
      rmse = std::sqrt(rmse / count);
      HERMES_LOG_VARIABLE(rmse);
      HERMES_LOG_VARIABLE(acc / count);
      csv.appendLine(methodName(args.interpolation_method), ", ", angle, ", ",
                     cell_size, ", ", rmse);
    }
    HERMES_LOG_VARIABLE((args.output_path / "data.csv").fullName());
    hermes::FileSystem::appendToFile(args.output_path / "data.csv", csv.str());
  }
  std::cerr << hermes::profiler::Profiler::dump();
  return 0;
}

int run() {
  SimPath sim_path(args.path);

  struct FrameData {
    hermes::Path frame;
    hermes::Path output_frame_path;
  };

  size_t sequential_frame_number = 1;
  size_t frame_number = 0;
  std::vector<FrameData> frame_data;
  for (auto &frame : sim_path.frames()) {
    if (frame.name() == "0")
      continue;
    frame_number++;

    if (frame_number % args.every_n != 0)
      continue;

    if (args.single_frame > 0 && frame_number != args.single_frame)
      continue;

    auto frame_time = std::stof(frame.name());
    auto frame_name =
        hermes::Str::concat(std::setfill('0'), std::setw(6), frame_time * 100);
    if (args.renumber_frames) {
      // use the dt to compute the actual frame number
      auto fn = hermes::Str::concat(std::setfill('0'), std::setw(6),
                                    frame_time / args.dt);
      HERMES_LOG("renaming frame {} to {}", frame_name, fn);
      frame_name = fn;
    } else if (args.renumber_frames_seq) {
      auto fn = hermes::Str::concat(std::setfill('0'), std::setw(6),
                                    sequential_frame_number);
      HERMES_LOG("renaming frame {} to {}", sequential_frame_number, fn);
      frame_name = fn;
    }

    auto output_frame_path = args.output_path / (frame_name + ".vdb");
    if (!args.force_write && args.single_frame != frame_number &&
        hermes::FileSystem::fileExists(output_frame_path))
      continue;

    frame_data.push_back(
        {.frame = frame, .output_frame_path = output_frame_path});

    sequential_frame_number++;
  }

  HERMES_LOG("Retrieving mesh ... ");
  PolyMesh mesh;
  mesh.load(args.path / "constant/polyMesh");
  HERMES_LOG("completed");

  // setup vdb
  VDBVolumeInterface *volume = nullptr;
  switch (args.interpolation_method) {
  case InterpolationMethod::HRBF:
    volume = new PUVolume<HRBFVolume>(mesh);
    break;
  case InterpolationMethod::CELL:
    volume = new PUVolume<CellVolume>(mesh);
    break;
  default:
    return -1;
  }

  // VecPUVolume u_volume(mesh);

  for (size_t i = 0; i < frame_data.size(); ++i) {
    const auto &frame = frame_data[i].frame;
    const auto &output_frame_path = frame_data[i].output_frame_path;
    // read frame data
    OpenFoamDict alpha_dict(frame / "alpha.snow");
    auto values = OpenFoamDict::parseValuesFrom<double>(
        alpha_dict["internalField"].value);
    if (values.size() != mesh.cells.size())
      values.resize(mesh.cells.size(), 0);
    // read velocity
    OpenFoamDict u_dict(frame / "U");
    auto u_values = OpenFoamDict::parseValuesFrom<hermes::vec3>(
        u_dict["internalField"].value);
    if (u_values.size() != mesh.cells.size())
      u_values.resize(mesh.cells.size(), hermes::vec3());
    // read grad
    OpenFoamDict g_dict(frame / "gradAlpha1");
    auto g_values = OpenFoamDict::parseValuesFrom<hermes::vec3>(
        g_dict["internalField"].value);
    if (g_values.size() != mesh.cells.size())
      g_values.resize(mesh.cells.size(), hermes::vec3());
    HERMES_LOG("converting {} with {} cells", frame, values.size());
    volume->update(i, values, g_values);
    // u_volume.update(i, u_values);
    auto vdb = volume->output(i, output_frame_path);
  }
  return 0;
}

int main(int argc, const char **argv) {
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  openvdb::initialize();

  // PUVolume<HRBFVolume>::testBlend();
  //  return 0;

  if (!args.read(argc, argv))
    return -1;
  if (args.test)
    return test();
  /*
else if (args.test_section) {
  hermes::Str csv;
  csv.appendLine("method, put_size, rmse");
  // args.debug = true;
  // args.interpolation_method = InterpolationMethod::NEAREST;
  // csv.appendLine(methodName(args.interpolation_method), ", ",
  // args.ptu_size,
  //                ", ", test_section());
  //  args.interpolation_method = InterpolationMethod::SHEPARD;
  //  csv.appendLine(methodName(args.interpolation_method), ", ",
  //  args.ptu_size,
  //                ", ", test_section());
  args.interpolation_method = InterpolationMethod::RBF;
  csv.appendLine(methodName(args.interpolation_method), ", ", args.ptu_size,
                 ", ", test_section());
  // args.interpolation_method = InterpolationMethod::HRBF;
  // csv.appendLine(methodName(args.interpolation_method), ", ",
  // args.ptu_size,
  //                ", ", test_section());
  HERMES_LOG_VARIABLE((args.output_path / "data.csv").fullName());
  hermes::FileSystem::writeFile(args.output_path / "data.csv", csv.str());
  return 0;
}*/
  return run();
}
