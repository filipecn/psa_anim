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

#include <dsl_mesh.h>
#include <hermes/common/arg_parser.h>
#include <hermes/common/debug.h>
#include <hermes/common/file_system.h>
#include <hermes/common/index.h>
#include <hermes/geometry/transform.h>
#include <hermes/logging/logging.h>
#include <hermes/numeric/numeric.h>
#include <mesh_utils.h>
#include <openfoam_poly_mesh.h>
#include <sim_path.h>

class Stats {
public:
  Stats(PolyMesh &m, const hermes::Path &o, float dx, float min_cut)
      : min_cut{min_cut}, rdx{1 / dx}, mesh{m}, output_path{o} {
    m.computeCellCenters();
    topology = psa_anim::computeFaceTopology(mesh.faces, mesh.cells);
    auto terrain = mesh.patch("terrain");
    HERMES_ASSERT(terrain.size > 0);
    // build terrain_height
    for (size_t face_id = 0; face_id < terrain.size; ++face_id) {
      auto center = mesh.faceCenter(terrain.start + face_id);
      auto key = encodeX(center.x);
      terrain_height[key] = center.z;
      terrain_x[key] = center.x;
    }
  }

  void storeInjectionVelocity(const hermes::Path &frame_path,
                              const std::string &frame_name) {
    OpenFoamDict dict(frame_path / "Uinj");
    auto values = OpenFoamDict::parseValuesFrom<hermes::vec3>(
        dict["boundaryField"]["terrain"]["value"].value);
    auto patch = mesh.patch("terrain");
    hermes::Str s;
    for (size_t face_id = 0; face_id < values.size(); ++face_id) {
      auto center = mesh.faceCenter(patch.start + face_id);
      auto p = center + 10.0f * values[face_id];
      s.appendLine("v ", p.x, " ", p.y, " ", p.z);
    }
    s.append("l");
    for (size_t face_id = 0; face_id < values.size(); ++face_id)
      s.append(" ", face_id + 1);
    s.appendLine();
    auto output_coords = output_path / (frame_name + "_injection.obj");
    hermes::FileSystem::writeFile(output_coords, s.str());
  }

  void storePlumeProfile(const hermes::Path &frame_path,
                         const std::string &frame_name) {
    OpenFoamDict dict(frame_path / "alpha.snow");
    auto values =
        OpenFoamDict::parseValuesFrom<double>(dict["internalField"].value);
    if (values.size() != mesh.cells.size())
      values.resize(mesh.cells.size(), 0);
    std::unordered_map<int, float> plume_height;
    for (size_t cell_id = 0; cell_id < values.size(); cell_id++) {
      auto center = mesh.cell_centers[cell_id];
      auto key = encodeX(center.x);

      if (values[cell_id] < 0.001)
        continue;

      if (!plume_height.count(key)) {
        plume_height[key] = center.z;
      }
      HERMES_ASSERT(terrain_height.count(key));
      plume_height[key] = std::max(terrain_height[key], center.z);
    }

    // build graph
    {
      std::vector<int> keys;
      for (const auto &item : plume_height)
        keys.emplace_back(item.first);
      std::sort(keys.begin(), keys.end());

      auto output_coords = output_path / (frame_name + "_plumes.obj");
      hermes::Str s;
      for (auto key : keys) {
        s.appendLine("v ", terrain_x[key], " 0.0 ", plume_height[key]);
      }
      s.append("l");
      for (size_t face_id = 0; face_id < keys.size(); ++face_id)
        s.append(" ", face_id + 1);
      s.appendLine();
      hermes::FileSystem::writeFile(output_coords, s.str());
    }
    {
      std::vector<int> keys;
      for (const auto &item : terrain_height)
        keys.emplace_back(item.first);
      std::sort(keys.begin(), keys.end());
      auto output_coords = output_path / (frame_name + "_plumes.csv");
      hermes::Str s;
      for (auto key : keys) {
        s.appendLine(terrain_x[key], " ",
                     std::max(0.f, plume_height[key] - terrain_height[key]));
      }
      hermes::FileSystem::writeFile(output_coords, s.str());
    }
  }

  void storeFront(const hermes::Path &frame_path,
                  const std::string &frame_name) {
    std::vector<double> dsl_h;
    {
      OpenFoamDict dict(frame_path / "dsl_h");
      dsl_h = OpenFoamDict::parseValuesFrom<double>(
          dict["boundaryField"]["terrain"]["value"].value);
    }
    std::vector<hermes::vec3> dsl_u;
    {
      OpenFoamDict dict(frame_path / "dsl_U");
      dsl_u = OpenFoamDict::parseValuesFrom<hermes::vec3>(
          dict["boundaryField"]["terrain"]["value"].value);
    }
    OpenFoamDict dict(frame_path / "dsl_front_sdf");
    auto values = OpenFoamDict::parseValuesFrom<double>(
        dict["boundaryField"]["terrain"]["value"].value);
    auto patch = mesh.patch("terrain");
    struct XAndInd {
      float x;
      size_t i;
    };
    std::vector<XAndInd> sorted_cells;
    float velocity = 0;
    float max_x = -1.f;
    for (size_t face_id = 0; face_id < values.size(); ++face_id) {
      auto center = mesh.faceCenter(patch.start + face_id);
      // auto v = dsl_u[face_id].length();
      // if (v > velocity) {
      //   velocity = v;
      //   max_x = center.x;
      // }
      if (dsl_h[face_id] > 0.1 && values[face_id] > 0.1 &&
          values[face_id] < 10000) {
        sorted_cells.push_back({center.x, face_id});
      }
    }
    if (!sorted_cells.empty())
      std::sort(
          sorted_cells.begin(), sorted_cells.end(),
          [](const XAndInd &a, const XAndInd &b) -> bool { return a.x > b.x; });
    if (!sorted_cells.empty()) {
      max_x = sorted_cells.front().x;
      int den = 0;
      for (int i = 0; i < 5 && i < sorted_cells.size(); ++i) {
        den++;
        velocity += dsl_u[sorted_cells[i].i].length();
      }
      velocity /= den;
    }
    {
      auto output_coords = output_path / (frame_name + "_front.obj");
      hermes::Str s;
      s.appendLine("v ", max_x, " 0.0 ", terrain_height[encodeX(max_x)]);
      hermes::FileSystem::writeFile(output_coords, s.str());
    }
    front_data.push_back({max_x, velocity});
  }

  struct FrontData {
    float x;
    float v;
  };
  std::vector<FrontData> front_data;

private:
  int encodeX(float x) { return x * rdx; }

  float min_cut{0.001};
  float rdx{0.5f};
  const PolyMesh &mesh;
  const hermes::Path &output_path;
  std::unordered_map<size_t, std::pair<i64, i64>> topology;
  std::unordered_map<int, float> terrain_height;
  std::unordered_map<int, float> terrain_x;
};

int main(int argc, const char **argv) {
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  hermes::ArgParser parser("stats", "profiles a given psl simulation");
  parser.addArgument("-i", "OpenFOAM simulation path", true);
  parser.addArgument("-o", "output path", true);
  parser.addArgument("--renumber-frames",
                     "Re-number frames based on the given dt");
  parser.addArgument("--renumber-frames-sequentially",
                     "Re-number frames sequentially");
  parser.addArgument("--dt",
                     "Re-number frames based on a given dt [default = 0.05]");
  parser.addArgument("--dx", "Cell size in x-axis");
  parser.addArgument("--min", "Filter values bellow this [default = 0.001]");
  parser.parse(argc, argv, true);

  hermes::Path path(parser.get<std::string>("-i"));
  hermes::Path output_path(parser.get<std::string>("-o"));
  bool renumber_frames = parser.check("--renumber-frames");
  bool renumber_frames_seq = parser.check("--renumber-frames-sequentially");
  auto dt = parser.get<float>("--dt", 0.05);
  auto dx = parser.get<float>("--dx", 2);
  const auto min_cut = parser.get<float>("--min", 0.001);

  // check output directory
  if (!output_path.isDirectory()) {
    HERMES_LOG_ERROR("invalid output path");
    return -1;
  }

  PolyMesh mesh;
  mesh.load(path / "constant/polyMesh");
  Stats stats(mesh, output_path, dx, min_cut);

  SimPath sim_path(path);

  size_t sequential_frame_number = 1;
  size_t frame_number = 0;

  std::vector<float> frame_times;
  for (auto &frame : sim_path.frames()) {
    if (frame.name() == "0")
      continue;
    frame_number++;

    auto frame_time = std::stof(frame.name());
    frame_times.emplace_back(frame_time);
    auto frame_name =
        hermes::Str::concat(std::setfill('0'), std::setw(6), frame_time * 100);

    if (renumber_frames) {
      // use the dt to compute the actual frame number
      auto fn =
          hermes::Str::concat(std::setfill('0'), std::setw(6), frame_time / dt);
      frame_name = fn;
    } else if (renumber_frames_seq) {
      auto fn = hermes::Str::concat(std::setfill('0'), std::setw(6),
                                    sequential_frame_number);
      frame_name = fn;
    }

    stats.storePlumeProfile(frame, frame_name);
    stats.storeFront(frame, frame_name);
    stats.storeInjectionVelocity(frame, frame_name);

    sequential_frame_number++;
  }

  auto front_file = output_path / "front_data.csv";
  HERMES_LOG_ERROR("writting {}: {} - {}", front_file, frame_times.size(),
                   stats.front_data.size());
  hermes::Str s;
  for (size_t i = 0; i < frame_times.size(); ++i) {
    s.appendLine(frame_times[i], ",", stats.front_data[i].x, ",",
                 stats.front_data[i].v);
  }
  hermes::FileSystem::writeFile(front_file, s.str());

  return 0;
}
