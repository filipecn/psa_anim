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
///\file psa_animpy.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-03-23
///
///\brief

#include <block_mesh.h>
#include <colors.h>
#include <dsl_mesh.h>
#include <hermes/geometry/point.h>
#include <hermes/logging/logging.h>
#include <hermes/numeric/numeric.h>
#include <msh_reader.h>
#include <openfoam_parser.h>
#include <openfoam_poly_mesh.h>
#include <post_process.h>
#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <sim_path.h>

#include <iomanip>
#include <sstream>
#include <utility>

namespace py = pybind11;

// stdout redirection
void printCallback(const hermes::Str &m) {
  py::print(m.str());
  std::cout << std::flush;
}

void initHermes(bool verbose = false) {
  std::cout << std::unitbuf;
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  hermes::Log::info_callback = printCallback;
  if (verbose) {
    hermes::Log::log_callback = [](const hermes::Str &m,
                                   hermes::logging_options o) {
      py::print(m.str());
      std::cout << std::flush;
    };
  } else {
    hermes::Log::log_callback = [](const hermes::Str &m,
                                   hermes::logging_options o) {};
  }
}

struct PSLDebugData {
  std::vector<int> global_face_id{};
  std::vector<int> local_face_id{};
  std::vector<double> dsl_area{};
  std::vector<double> snow_cover_mass{};
  std::vector<double> dm{};
  std::vector<double> entrainment_mass{};
  std::vector<double> front_weight{};
  std::vector<double> dsl_velocity_mag{};
  std::vector<double> dsl_height{};
  std::vector<double> we{};
  std::vector<double> slope_angle{};
};

///
/// \param stl_path
/// \param obj_path
/// \param a
/// \param b
/// \param cell_size
/// \param width
void stl2quad(const std::string &stl_path, const std::string &obj_path,
              const py::array_t<float> &a, const py::array_t<float> &b,
              float cell_size, float width, bool twoDim, bool verbose = false) {
  initHermes(verbose);
  HERMES_LOG_VARIABLE(stl_path);
  HERMES_LOG_VARIABLE(obj_path);
  auto n = a.unchecked<1>();
  HERMES_ASSERT(n.size() == 2);
  hermes::point2 axis_a{n[0], n[1]};
  auto n2 = b.unchecked<1>();
  HERMES_ASSERT(n2.size() == 2);
  hermes::point2 axis_b{n2[0], n2[1]};
  hermes::vec2 axis = hermes::normalize(axis_b - axis_a);
  hermes::point2 origin = axis_a; // + axis * 100;
  auto terrain = psa_anim::loadSTL(stl_path, "terrain");
  auto grid = psa_anim::extractGrid(*terrain, axis, axis_a,
                                    {(axis_b - origin).length() * 0.95f, width},
                                    cell_size, twoDim);
  psa_anim::info(grid);
  psa_anim::saveOBJ(grid, obj_path);
}

///
/// \param input
/// \param output
/// \param deep_height
/// \param bottom_type 0 - XY plane, 1 - inclined, 2 - extrusion
/// \param close_top
void terrainBase(const std::string &input, const std::string &output,
                 float deep_height, int bottom_type, bool close_top,
                 bool verbose = false) {
  initHermes(verbose);
  HERMES_LOG("Creating terrain base parameters");
  HERMES_LOG("==================================================");
  HERMES_LOG_VARIABLE(input);
  HERMES_LOG_VARIABLE(output);
  HERMES_LOG_VARIABLE(deep_height);
  HERMES_LOG_VARIABLE(bottom_type);
  HERMES_LOG_VARIABLE((int)close_top);
  HERMES_LOG("==================================================");
  HERMES_ASSERT(bottom_type <= 2)
  auto terrain = psa_anim::loadOBJ(input);
  HERMES_ASSERT(terrain)
  auto edge_topo = psa_anim::computeEdgeTopology(terrain->faces);

  // copy data
  auto base = *terrain;

  // in order to create walls we will need the boundary edges
  std::vector<std::pair<size_t, size_t>> boundary_edges;
  for (const auto &e : edge_topo) {
    auto edge_key = e.first;
    auto cell_a = e.second.first;
    auto cell_b = e.second.second;
    HERMES_ASSERT(cell_a >= 0) // safe check on the input
    if (cell_b < 0)
      // edge is boundary
      boundary_edges.emplace_back(edge_key);
  }

  // if bottom mode is just a extrusion then we just duplicate everything
  // the bottom plane is defined by its normal and a point based on the minimal
  // z terrain coordinate here we get the deepest point in the terrain
  auto deepest_point = base.vertices.front();
  for (const auto &v : base.vertices)
    if (v.z < deepest_point.z)
      deepest_point = v;

  // if the bottom_type is zero, then the normal is on the -z direction
  if (bottom_type == 0 || bottom_type == 1) {
    hermes::normal3 normal(0, 0, -1);
    if (bottom_type == 1) {
      // if the bottom is inclined, the the normal is based on the mean normal
      // of the terrain
      normal = {0, 0, 0};
      for (const auto &face : terrain->faces)
        normal += psa_anim::faceNormal(terrain->vertices, face);
      normal /= terrain->faces.size();
      normal = hermes::normalize(normal);
    }
    hermes::Plane plane(-normal,
                        deepest_point - hermes::vec3(0, 0, deep_height));
    // project vertices to bottom plane
    for (size_t i = 0; i < terrain->vertices.size(); ++i) {
      hermes::Ray3 r(base.vertices[i], hermes::vec3(0, 0, -1));
      auto t = hermes::GeometricPredicates::intersect(plane, r);
      HERMES_ASSERT(t)
      base.vertices[i] = r(*t);
    }
  } else if (bottom_type == 2) {
    // in the case of the extruded bottom we just translate the vertices
    for (size_t i = 0; i < terrain->vertices.size(); ++i)
      base.vertices[i].z -= deep_height;
  }

  // if the top is closed we duplicate all faces
  if (close_top) {
    HERMES_NOT_IMPLEMENTED
  } else {
    // duplicate boundary vertices
    //                 input index -> output index
    std::unordered_map<size_t, size_t> boundary_index_map;
    for (const auto &edge : boundary_edges) {
      if (!boundary_index_map.count(edge.first)) {
        boundary_index_map[edge.first] = base.vertices.size();
        base.vertices.emplace_back(terrain->vertices[edge.first]);
      }
      if (!boundary_index_map.count(edge.second)) {
        boundary_index_map[edge.second] = base.vertices.size();
        base.vertices.emplace_back(terrain->vertices[edge.second]);
      }
    }
    // flip bottom faces (here we only have the bottom, no walls yet
    //    for (auto &face : base.faces)
    //      psa_anim::flipOrder(face);
    // create walls
    for (const auto &edge : boundary_edges) {
      auto e = edge_topo.find(edge);
      HERMES_ASSERT(e != edge_topo.end())
      HERMES_ASSERT(e->second.second < 0)
      HERMES_ASSERT(e->second.first >= 0 && e->second.first < base.faces.size())

      auto a = edge.first;
      auto b = edge.second;
      HERMES_ASSERT(boundary_index_map.count(a))
      HERMES_ASSERT(boundary_index_map.count(b))
      auto c = boundary_index_map[b];
      auto d = boundary_index_map[a];

      if (base.edgeIsSorted(e->second.first, edge.first, edge.second)) {
        std::swap(a, b);
        std::swap(c, d);
      }
      base.faces.push_back({a, b, c, d});
    }
  }
  // output mesh
  HERMES_LOG("Storing base terrain in [{}]", output);
  psa_anim::saveOBJ(base, output);
}

/**********************************************************************************************************************/
/*                                                                                                BlockMeshTerrain_py
 */
/* Quad surface mesh to block mesh description converter. */
/**********************************************************************************************************************/
class BlockMeshTerrain_py {
public:
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  ///
  BlockMeshTerrain_py(bool verbose = false) { initHermes(verbose); }
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  ///
  /// \param patch_name
  /// \param normal
  void addPatchDirection(const std::string &patch_name,
                         const py::array_t<float> &normal) {
    auto n = normal.unchecked<1>();
    HERMES_ASSERT(n.size() == 3);
    boundaries_[patch_name].emplace_back(n[0], n[1], n[2]);
  }
  ///
  /// \param patch_name
  /// \param type
  void setBoundaryType(const std::string &patch_name, const std::string &type) {
    boundary_types_[patch_name] = type;
  }
  ///
  /// \param path
  void loadOBJ(const std::string &path) {
    HERMES_LOG_VARIABLE(path);
    auto mesh = psa_anim::loadOBJ(path);
    HERMES_ASSERT(mesh);
    HERMES_LOG_VARIABLE(mesh->faces.size());
    HERMES_LOG_VARIABLE(mesh->vertices.size());
    bm_ = psa_anim::blockMesh(mesh->vertices, mesh->faces, boundaries_, height_,
                              config_, slope_direction_, slope_point_,
                              top_is_inclined_);
    for (const auto &b : boundary_types_)
      bm_.setFaceType(b.first, b.second);
  }
  ///
  /// \param path
  void save(const std::string &path) {
    HERMES_LOG("storing blockMeshDict at {}", path);
    /*
        psa_anim::FaceMesh m;
        std::vector<size_t> vs = {
            412, 411, 0, 1, 50143, 50142, 49731, 49732
        };
        hermes::Str s;
        for (auto v : vs)
          s.appendLine("v ", bm_.vertexAt(v).x, " ",
                       bm_.vertexAt(v).y, " ",
                       bm_.vertexAt(v).z);
        std::vector<std::vector<size_t>> indices = {
            {0, 4, 7, 3}, // f0
            {1, 2, 6, 5}, // f1
            {0, 1, 5, 4}, // f2
            {2, 3, 7, 6}, // f3
            {0, 3, 2, 1}, // f4
            {4, 5, 6, 7}, // f5
        };

        hermes::vec3 a = bm_.vertexAt(vs[indices[4][1]]) -
            bm_.vertexAt(vs[indices[4][0]]);
        hermes::vec3 b = bm_.vertexAt(vs[indices[4][2]]) -
            bm_.vertexAt(vs[indices[4][1]]);
        HERMES_LOG_VARIABLE(hermes::normalize(hermes::cross(a,b)));
        for (auto &in : indices)
          for (auto &v : in)
            v += 1;
        for (int i = 0; i < 6; ++i)
          s.appendLine("f ", hermes::Str::join(indices[i], " "));

        hermes::FileSystem::writeFile("/home/filipecn/Desktop/debug.obj",
       s.str()); bm_.saveOBJ("/home/filipecn/Desktop/debug.obj");
    */

    bm_.save(path);
  }
  void setHeight(real_t height) { height_ = height; }
  void setSlopeDirection(const py::array_t<float> &v) {
    auto n = v.unchecked<1>();
    HERMES_ASSERT(n.size() == 2);
    slope_direction_.x = n[0];
    slope_direction_.y = n[1];
    slope_direction_ = hermes::normalize(slope_direction_);
  }
  void setSlopePoint(const py::array_t<float> &v) {
    auto n = v.unchecked<1>();
    HERMES_ASSERT(n.size() == 2);
    slope_point_.x = n[0];
    slope_point_.y = n[1];
  }
  void setTopIsInclined(bool b) { top_is_inclined_ = b; }
  void setFaceMatching(bool b) { config_.face_matching = b; }
  void setGrading(const py::array_t<float> &v) {
    auto n = v.unchecked<1>();
    HERMES_ASSERT(n.size() == 3);
    config_.grading.x = n[0];
    config_.grading.y = n[1];
    config_.grading.z = n[2];
  }
  void setBlockResolution(const py::array_t<int> &v) {
    auto n = v.unchecked<1>();
    HERMES_ASSERT(n.size() == 3);
    config_.resolution.width = n[0];
    config_.resolution.height = n[1];
    config_.resolution.depth = n[2];
  }

private:
  std::map<std::string, std::vector<hermes::vec3>> boundaries_;
  std::map<std::string, std::string> boundary_types_;
  psa_anim::BlockMeshDescription bm_;
  real_t height_{50};
  hermes::vec2 slope_direction_{1, 0};
  hermes::point2 slope_point_{0, 0};
  bool top_is_inclined_{true};
  psa_anim::BlockMeshConfig config_;
};

/**********************************************************************************************************************/
/*                                                                                                             DSL_py
 */
/* Set of algorithms to process the DSL simulation output for PSL simulation
 * input.                                   */
/* The following data is generated for each simulation step: */
/*  - front sdf: front distance field defined in every surface cell */
/*  - slope: terrain slope computed in every surface cell */
/*  - dsl height field */
/*  - dsl velocity field */
/*  - front cells */
/*                                                                                                                    */
/* Output data is produced for surface face center points, thus output geometry
 * is a set of points (and not a mesh).  */
/**********************************************************************************************************************/
class DSL_py {
public:
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  ///
  DSL_py(bool verbose = false) { initHermes(verbose); }
  ///
  /// \param poly_mesh openfoam mesh object
  /// \param patch_id openfoam mesh patch identified
  /// \param sim_path directory containing the DSL simulation steps
  DSL_py(const PolyMesh *poly_mesh, const std::string &patch_name,
         const SimPath &sim_path, bool verbose = false)
      : sim_path_(sim_path) {
    initHermes(verbose);
    dsl_mesh_.set(poly_mesh, patch_name);
  }
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  //                                                                                                           access
  ///
  /// \return
  [[nodiscard]] std::string patchName() const { return dsl_mesh_.patchName(); }
  ///
  /// \return
  [[nodiscard]] size_t size() const { return dsl_mesh_.size(); }
  // *******************************************************************************************************************
  //                                                                                                               IO
  // *******************************************************************************************************************
  // dsl specifics
  ///
  /// \param sim_path
  /// \param output_path
  /// \param jitter_points
  void preparePSLInput(std::string output_path, bool jitter_points = false,
                       bool empty_dsl = false) {
    HERMES_LOG("Preparing psl input");
    HERMES_LOG_AND_RETURN_IF_NOT(dsl_mesh_.good(), "invalid patch id");
    HERMES_LOG_VARIABLE(output_path);
    hermes::Path path(std::move(output_path));
    writePointsFile(path / "points", jitter_points);

    sim_path_.setCurrentFrame(1);

    for (; !sim_path_.finished(); ++sim_path_) {
      // load frame data
      dsl_mesh_.load(sim_path_);
      auto frame = sim_path_.currentFramePath();
      auto output_frame_path = path / frame.name();
      if (!hermes::FileSystem::isDirectory(output_frame_path)) {
        HERMES_ASSERT(hermes::FileSystem::mkdir(output_frame_path));
      }

      // TODO hack: here we synthesize the frame 0 by copying the second frame
      // with zero velocity
      if (sim_path_.currentFrame() == 1) {
        auto zero_frame_path = path / "0";
        if (!hermes::FileSystem::isDirectory(zero_frame_path)) {
          HERMES_ASSERT(hermes::FileSystem::mkdir(zero_frame_path));
        }
        extractField("volScalarField", zero_frame_path / "dsl_W", "fa_W");
        extractField("volScalarField", zero_frame_path / "dsl_divU",
                     "fa_divUs");
        extractField("volScalarField", zero_frame_path / "dsl_h", "fa_h");
        extractField("volVectorField", zero_frame_path / "dsl_U", "fa_Us",
                     true);
        writeFrontDistanceField(zero_frame_path / "dsl_front_sdf",
                                zero_frame_path / "front");
        computeSlopeField(zero_frame_path / "dsl_slope");
      }

      extractField("volScalarField", output_frame_path / "dsl_W", "fa_W",
                   empty_dsl);
      extractField("volScalarField", output_frame_path / "dsl_divU", "fa_divUs",
                   empty_dsl);
      extractField("volScalarField", output_frame_path / "dsl_h", "fa_h",
                   empty_dsl);
      extractField("volVectorField", output_frame_path / "dsl_U", "fa_Us",
                   empty_dsl);
      writeFrontDistanceField(output_frame_path / "dsl_front_sdf",
                              output_frame_path / "front", empty_dsl);
      computeSlopeField(output_frame_path / "dsl_slope");
    }
  }
  ///
  /// \param output_path
  /// \param max_height
  /// \param with_height_map
  /// \param renumber_frames
  /// \param only_missing_frames
  /// \param force_write
  /// \param single_frame_number
  void extractSurface(std::string output_path, float max_height = 10,
                      bool with_height_map = false,
                      bool with_distance_map = false,
                      bool renumber_frames = false,
                      bool only_missing_frames = false,
                      bool force_write = false, int single_frame_number = -1) {
    HERMES_LOG_AND_RETURN_IF_NOT(dsl_mesh_.good(), "invalid patch id");
    HERMES_LOG("Extracting DSL surface");

    float dt = 0.05;

    hermes::Path path(std::move(output_path));

    auto palette = psa_anim::ColorPalettes::Batlow();

    sim_path_.setCurrentFrame(1);
    size_t output_count = 0;
    for (; !sim_path_.finished(); ++sim_path_) {
      auto frame = sim_path_.currentFramePath();
      // output
      auto frame_time = std::stof(frame.name());
      auto frame_name = hermes::Str::concat(std::setfill('0'), std::setw(6),
                                            frame_time * 100);

      if (renumber_frames)
        frame_name = hermes::Str::concat(std::setfill('0'), std::setw(6),
                                         frame_time / dt);
      auto frame_output_path = path / (frame_name + ".ply");

      if (single_frame_number > 0 &&
          single_frame_number != sim_path_.currentFrame())
        continue;

      if (frame_output_path.exists() && !force_write && single_frame_number < 0)
        continue;

      if (frame_output_path.exists() && only_missing_frames)
        continue;

      // load frame data
      dsl_mesh_.load(sim_path_);
      std::vector<double> front_distance_values;
      std::vector<double> ws;
      std::vector<double> mag_us;
      std::vector<double> u_inj;
      auto dsl_surface_mesh =
          dsl_mesh_.surfaceMesh(front_distance_values, ws, mag_us, u_inj);
      float max_distance = 3000;

      // compute vertices colors
      std::vector<psa_anim::Color> colors(dsl_surface_mesh.vertices.size());
      std::vector<psa_anim::Color> sdf_colors(dsl_surface_mesh.vertices.size());
      std::vector<float> height_values(dsl_surface_mesh.vertices.size());
      auto bottom_vertex_count = dsl_surface_mesh.vertices.size() / 2;
      for (size_t i = 0; i < bottom_vertex_count; ++i) {
        height_values[i] = height_values[i + bottom_vertex_count] =
            hermes::distance(
                dsl_surface_mesh.vertices[i],
                dsl_surface_mesh.vertices[i + bottom_vertex_count]);
        colors[i] = colors[i + bottom_vertex_count] =
            palette(hermes::Numbers::clamp(
                        hermes::distance(
                            dsl_surface_mesh.vertices[i],
                            dsl_surface_mesh.vertices[i + bottom_vertex_count]),
                        0.f, max_height) /
                    max_height);
        sdf_colors[i] = sdf_colors[i + bottom_vertex_count] =
            palette(hermes::Numbers::clamp((float)front_distance_values[i], 0.f,
                                           max_distance) /
                    max_distance);
      }

      // store in ply format!
      hermes::Str content;
      content.appendLine("ply");
      content.appendLine("format ascii 1.0");
      content.appendLine("element vertex ", dsl_surface_mesh.vertices.size());
      content.appendLine("property float x");
      content.appendLine("property float y");
      content.appendLine("property float z");
      if (with_height_map) {
        content.appendLine("property uchar red");
        content.appendLine("property uchar green");
        content.appendLine("property uchar blue");
      }
      if (with_distance_map) {
        content.appendLine("property uchar distred");
        content.appendLine("property uchar distgreen");
        content.appendLine("property uchar distblue");
      }
      content.appendLine("property float h");
      content.appendLine("property float front");
      content.appendLine("property float u");
      content.appendLine("property float w");
      content.appendLine("property float inj");
      content.appendLine("element face ", dsl_surface_mesh.faces.size());
      content.appendLine("property list uchar int vertex_index");
      content.appendLine("end_header");
      for (size_t i = 0; i < dsl_surface_mesh.vertices.size(); ++i) {
        content.append(dsl_surface_mesh.vertices[i].x, " ",
                       dsl_surface_mesh.vertices[i].y, " ",
                       dsl_surface_mesh.vertices[i].z, " ");
        if (with_height_map) {
          content.append(static_cast<size_t>(colors[i].r * 255.0), " ",
                         static_cast<size_t>(colors[i].g * 255.0), " ",
                         static_cast<size_t>(colors[i].b * 255.0));
        }
        if (with_distance_map) {
          content.append(static_cast<size_t>(sdf_colors[i].r * 255.0), " ",
                         static_cast<size_t>(sdf_colors[i].g * 255.0), " ",
                         static_cast<size_t>(sdf_colors[i].b * 255.0));
        }
        content.append(" ", height_values[i % bottom_vertex_count]);
        content.append(" ", front_distance_values[i % bottom_vertex_count]);
        content.append(" ", mag_us[i % bottom_vertex_count]);
        content.append(" ", ws[i % bottom_vertex_count]);
        content.append(" ", u_inj[i % bottom_vertex_count]);
        content.appendLine();
      }
      for (const auto &face : dsl_surface_mesh.faces)
        content.appendLine(face.size(), " ", hermes::Str::join(face, " "));

      if (output_count % 100 == 0)
        HERMES_LOG("storing dsl surface file [{}]", frame_output_path);
      hermes::FileSystem::writeFile(frame_output_path, content.str());
      output_count++;
    }
  }

protected:
  /// Writes OpenFOAM's points file containing the set of face centers
  /// \param path output file
  /// \param jitter_points options to jitter points in y (side) direction to
  /// avoid all collinear in a 2d example
  void writePointsFile(hermes::Path path, bool jitter_points) {
    hermes::Str str("/*--------------------------------*- C++ "
                    "-*----------------------------------*\n"
                    "| =========                 |                             "
                    "                    |\n"
                    "| \\\\      /  F ield         | OpenFOAM: The Open Source "
                    "CFD Toolbox           |\n"
                    "|  \\\\    /   O peration     | Version:  2012            "
                    "                      |\n"
                    "|   \\\\  /    A nd           | Website:  "
                    "www.openfoam.com                      |\n"
                    "|    \\\\/     M anipulation  |                           "
                    "                      |\n"
                    "\\*-------------------------------------------------------"
                    "--------------------*/\n"
                    "FoamFile\n"
                    "{\n"
                    "    version     2.0;\n"
                    "    format      ascii;\n"
                    "    class       pointField;\n"
                    "    location    \"constant/boundaryData/points\";\n"
                    "    object      points;\n"
                    "}\n"
                    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * "
                    "* * * * * * * * * * //\n\n\n");
    str.appendLine(dsl_mesh_.size());
    str.appendLine("(");
    // store points file (face centers)
    auto centers = dsl_mesh_.faceCenters();
    for (size_t i = 0; i < centers.size(); ++i) {
      auto dy = i % 2 ? -0.05 : 0.05;
      if (!jitter_points)
        dy = 0;
      str.appendLine(hermes::Str::format("({} {} {})", centers[i].x,
                                         centers[i].y + dy, centers[i].z));
    }
    str.appendLine(")\n");
    HERMES_LOG("writing points file: {}", path.fullName());
    hermes::FileSystem::writeFile(path, str.str());
  }
  /// Extract a field for the selected poly mesh patch and stores in a separate
  /// OpenFOAM file \param input file field input path \param field_type
  /// [volScalarField | volVectorField] \param output_path output file path
  /// \param field_name
  void extractField(std::string field_type, hermes::Path output_path,
                    std::string field_name, bool is_zero = false) {
    hermes::Str str(hermes::Str::format("FoamFile\n{\n"
                                        "    version     2.0;\n"
                                        "    format      ascii;\n"
                                        "    class       {};\n"
                                        "    location    \"{}\";\n"
                                        "    object       {};\n}\n\n",
                                        field_type, output_path.fullName(),
                                        field_name));
    auto n = dsl_mesh_.size();
    str.appendLine(n, "\n(");
    if (is_zero) {
      // assumes zero
      if (field_type == "volScalarField")
        for (size_t i = 0; i < n; ++i)
          str.appendLine(0);
      else if (field_type == "volVectorField")
        for (size_t i = 0; i < n; ++i)
          str.appendLine("(0 0 0)");
    } else {
      // read field
      if (field_type == "volScalarField") {
        for (auto value : dsl_mesh_.scalarField(field_name))
          str.appendLine(value);
      } else if (field_type == "volVectorField") {
        for (auto value : dsl_mesh_.vectorField(field_name))
          str.appendLine(
              hermes::Str::format("({} {} {})", value.x, value.y, value.z));
      }
    }
    str.appendLine(')');
    hermes::FileSystem::writeFile(output_path, str.str());
  }
  // *******************************************************************************************************************
  //                                                                                                        FRONT SDF
  // *******************************************************************************************************************
  /// Write front distance scalar field into an OpenFOAM dict file for a given
  /// input \param output_path output file path \param input_h DSL height field
  /// file \param input_U DSL velocity field file
  void writeFrontDistanceField(hermes::Path output_path,
                               hermes::Path output_front_path,
                               bool empty_dsl = false) {
    const auto &front = dsl_mesh_.frontCells();
    std::vector<double> dsl_front_sdf = dsl_mesh_.scalarField("front_sdf");
    // write front
    {
      hermes::Str str;
      for (const auto &front_cell : front)
        str.appendLine(
            front_cell.cell_id, " ",
            hermes::Str::join(std::vector<real_t>(
                                  {front_cell.position.x, front_cell.position.y,
                                   front_cell.position.z, front_cell.velocity.x,
                                   front_cell.velocity.y, front_cell.velocity.z,
                                   front_cell.height}),
                              " "));
      hermes::FileSystem::writeFile(output_front_path, str.str());
    }
    // write file
    hermes::Str str(hermes::Str::format("FoamFile\n{\n"
                                        "    version     2.0;\n"
                                        "    format      ascii;\n"
                                        "    class       volScalarField;\n"
                                        "    location    \"{}\";\n"
                                        "    object       dsl_slope;\n}\n\n",
                                        output_path.fullName()));
    str.appendLine(dsl_mesh_.size(), "\n(");
    for (auto sd : dsl_front_sdf)
      if (empty_dsl)
        str.appendLine(hermes::Numbers::greatest_f32());
      else
        str.appendLine(sd);
    str.appendLine(')');
    hermes::FileSystem::writeFile(output_path, str.str());
  }
  // *******************************************************************************************************************
  //                                                                                                            SLOPE
  // *******************************************************************************************************************
  /// Computes the terrain slope for each cell in the selected poly mesh patch
  /// \param output_path
  void computeSlopeField(hermes::Path output_path) {
    hermes::Str str(hermes::Str::format("FoamFile\n{\n"
                                        "    version     2.0;\n"
                                        "    format      ascii;\n"
                                        "    class       volScalarField;\n"
                                        "    location    \"{}\";\n"
                                        "    object       dsl_slope;\n}\n\n",
                                        output_path.fullName()));
    const auto &slope_angles = dsl_mesh_.scalarField("slope");
    str.appendLine(dsl_mesh_.size(), "\n(");
    for (auto angle : slope_angles)
      str.appendLine(angle);
    str.appendLine(')');
    hermes::FileSystem::writeFile(output_path, str.str());
  }

  DSLMesh dsl_mesh_;
  SimPath sim_path_;
};

/**********************************************************************************************************************/
/*                                                                                                   ProceduralDSL_py
 */
/* Set of procedural DSL data */
/**********************************************************************************************************************/
class ProceduralDSL_py : public DSL_py {
public:
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  ///
  ///
  ProceduralDSL_py(bool verbose = false) { initHermes(verbose); }
  ///
  /// \param poly_mesh openfoam mesh object
  /// \param patch_id openfoam mesh patch identified
  /// \param sim_path directory containing the DSL simulation steps
  ProceduralDSL_py(const PolyMesh *poly_mesh, const std::string &patch_name,
                   const SimPath &sim_path, bool verbose = false)
      : DSL_py(poly_mesh, patch_name, sim_path, verbose) {}
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  //                                                                                                           access
  ///
  /// \return
  [[nodiscard]] std::string patchName() const { return dsl_mesh_.patchName(); }
  ///
  /// \return
  [[nodiscard]] size_t size() const { return dsl_mesh_.size(); }
  //                                                                                                           config
  void setDuration(float duration) { duration_ = duration; }
  void setWriteInterval(float dt) { dt_ = dt; }
  void addRotatingDisk(float path_radius, float disk_radius, float velocity,
                       float cx, float cy) {
    rotating_disks_.push_back({.path_radius = path_radius,
                               .radius = disk_radius,
                               .velocity = velocity,
                               .path_center = {cx, cy, 0}});
  }
  void addMovingBar(float cx, float cy, float sx, float sy, float velocity) {
    moving_bars_.push_back({.center = {cx, cy, 0},
                            .half_size = {sx / 2, sy / 2, 10},
                            .velocity = velocity});
  }
  void addSteadyFront(float cx, float r, float v) {
    steady_fronts_.push_back({.radius = r, .center_x = cx, .velocity_x = v});
  }

  void write(const std::string &spath) {
    HERMES_LOG("Preparing procedural DSL output (PSL input)");
    HERMES_LOG_AND_RETURN_IF_NOT(dsl_mesh_.good(), "invalid patch id");
    hermes::Path path(std::move(spath));

    path /= dsl_mesh_.patchName();

    if (!hermes::FileSystem::isDirectory(path)) {
      HERMES_ASSERT(hermes::FileSystem::mkdir(path));
    }
    writePointsFile(path / "points", true);

    std::vector<double> fa_h(dsl_mesh_.size());
    std::vector<double> fa_W(dsl_mesh_.size());
    std::vector<double> fa_divUs(dsl_mesh_.size());
    std::vector<hermes::vec3> fa_Us(dsl_mesh_.size());

    size_t precision = 0;
    float dt = dt_;
    while (dt != 0 && int(dt / 10) == 0) {
      dt *= 10;
      precision++;
    }

    auto centers = dsl_mesh_.faceCenters();
    for (float time = 0; time <= duration_; time += dt_) {
      std::stringstream s;
      s << std::fixed << std::setprecision(precision) << time;
      auto frame_path = path / s.str();

      if (!hermes::FileSystem::isDirectory(frame_path)) {
        HERMES_ASSERT(hermes::FileSystem::mkdir(frame_path));
      }

      // setup height
      for (size_t i = 0; i < centers.size(); ++i) {
        for (const auto &disk : rotating_disks_)
          disk.eval(time, centers[i], fa_h[i], fa_Us[i], fa_W[i], fa_divUs[i]);
        for (const auto &bar : moving_bars_)
          bar.eval(time, centers[i], fa_h[i], fa_Us[i], fa_W[i], fa_divUs[i]);
        for (const auto &front : steady_fronts_)
          front.eval(time, centers[i], fa_h[i], fa_Us[i], fa_W[i], fa_divUs[i]);
      }

      dsl_mesh_.load(fa_h, fa_Us, fa_W, fa_divUs);

      extractField("volScalarField", frame_path / "dsl_h", "fa_h", false);
      extractField("volScalarField", frame_path / "dsl_W", "fa_W", false);
      extractField("volScalarField", frame_path / "dsl_divU", "fa_divUs",
                   false);
      extractField("volVectorField", frame_path / "dsl_U", "fa_Us", false);
      writeFrontDistanceField(frame_path / "dsl_front_sdf",
                              frame_path / "front", false);
      computeSlopeField(frame_path / "dsl_slope");
    }
  }

private:
  struct SteadyFront {
    float radius;
    float center_x;
    float velocity_x;
    void eval(float time, const hermes::point3 &p, double &h, hermes::vec3 &u,
              double &w, double &divU) const {
      h = 0;
      u = {0, 0, 0};
      divU = 0;
      w = 0;
      if (p.x < center_x ||
          hermes::distance(p, hermes::point3(center_x, 0, p.z)) > radius)
        return;
      h = 1;
      w = 1;
      u = {velocity_x, 0, 0};
    }
  };
  struct RotatingDisk {
    float path_radius;
    float radius;
    float velocity;
    hermes::point3 path_center;
    void eval(float time, const hermes::point3 &p, double &h, hermes::vec3 &u,
              double &w, double &divU) const {
      h = 0;
      u = {0, 0, 0};
      w = 0;
      divU = 0;
      auto angle = (velocity * time) / path_radius;
      hermes::point3 center(path_radius * std::cos(angle),
                            path_radius * std::sin(angle), 0);
      auto r = center - path_center;
      if (hermes::distance2(center + hermes::vec3(path_center), p) >
          radius * radius)
        return;
      h = 1;
      hermes::vec3 v(-center.y, center.x, 0);
      v.normalize();
      u = velocity * v;
    }
  };
  struct MovingBar {
    hermes::point3 center;
    hermes::vec3 half_size;
    float velocity;
    void eval(float time, const hermes::point3 &p, double &h, hermes::vec3 &u,
              double &w, double &divU) const {
      h = 0;
      u = {0, 0, 0};
      w = 0;
      divU = 0;

      auto cur_center = center;
      cur_center.x += velocity * time;

      hermes::BBox3 bbox(cur_center - half_size, cur_center + half_size);

      if (bbox.contains(p)) {
        u.x = velocity;
        h = 1;
      }
    }
  };
  float duration_;
  float dt_;
  std::vector<RotatingDisk> rotating_disks_;
  std::vector<MovingBar> moving_bars_;
  std::vector<SteadyFront> steady_fronts_;
};

/**********************************************************************************************************************/
/*                                                                                                              OF_py
 */
/* OpenFOAM's simulation data (including geometry and scalar/vector fields) */
/**********************************************************************************************************************/
class OF_py {
public:
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  ///
  OF_py(bool verbose = false) { initHermes(verbose); }
  /// Sets the simulation directory containing the mesh geometry in a
  /// constant/polyMesh sub-directory \param sim_path simulation directory path
  /// \return true if data is successfully loaded
  bool setSimPath(const std::string &sim_path) {
    HERMES_LOG("Load constant/polyMesh");
    HERMES_LOG("path: {}", sim_path);
    hermes::Path path(sim_path);
    mesh_.load(path / "constant/polyMesh");
    HERMES_LOG("Compute cell bounds");
    mesh_.computeCellBounds();
    HERMES_LOG("Compute cell centers");
    mesh_.computeCellCenters();
    HERMES_LOG("List fields");
    sim_path_ = SimPath(path);
    sim_path_.listFields(scalar_fields_, vector_fields_);
    return true;
  }
  /// Sets the directory containing the list of directories of the simulation
  /// steps \param frames_path \return true if data is successfully loaded
  bool setFramesPath(const std::string &frames_path) {
    sim_path_ = SimPath(frames_path);
    sim_path_.listFields(scalar_fields_, vector_fields_);
    return true;
  }
  // *******************************************************************************************************************
  //                                                                                                       SIMULATION
  // *******************************************************************************************************************
  ///
  /// \param frame_id
  void setCurrentFrame(size_t frame_id) {
    sim_path_.setCurrentFrame(frame_id);
    sim_path_.listFields(scalar_fields_, vector_fields_);
  }
  ///
  /// \return
  [[nodiscard]] std::string currentFramePath() const {
    return sim_path_.currentFramePath().fullName();
  }
  ///
  /// \return
  [[nodiscard]] size_t frameCount() const { return sim_path_.frames().size(); }
  // *******************************************************************************************************************
  //                                                                                                      DATA FIELDS
  // *******************************************************************************************************************
  ///
  /// \return
  [[nodiscard]] std::vector<std::string> scalarFields() const {
    return scalar_fields_;
  }
  /// Retrieves scalar field values for a given patch
  /// \param field_name
  /// \param patch_name
  /// \return list of values
  [[nodiscard]] std::vector<float>
  scalarFieldData(const std::string &field_name,
                  const std::string &patch_name) const {
    OpenFoamDict dict(sim_path_.currentFramePath() / field_name);
    // dict from patch
    const auto &boundary_patches = dict["boundaryField"];
    std::vector<float> values;
    auto patch_dict = boundary_patches.fields.find(patch_name);
    if (patch_dict != boundary_patches.fields.end()) {
      auto patch_values_dict = patch_dict->second.fields.find("value");
      if (patch_values_dict != patch_dict->second.fields.end()) {
        values =
            OpenFoamDict::parseValuesFrom<f32>(patch_values_dict->second.value);
      }
    }
    return values;
  }
  // *******************************************************************************************************************
  //                                                                                                          PATCHES
  // *******************************************************************************************************************
  [[nodiscard]] ProceduralDSL_py proceduralDSL(const std::string &patch_name,
                                               bool verbose = false) const {
    return ProceduralDSL_py(&mesh_, patch_name, sim_path_, verbose);
  }
  [[nodiscard]] DSL_py DSL(const std::string &patch_name,
                           bool verbose = false) const {
    return DSL_py(&mesh_, patch_name, sim_path_, verbose);
  }
  /// Exports patch geometry to OBJ format
  /// \param patch_name patch identifier
  /// \param output_path obj output file path
  void exportOBJ(const std::string &patch_name,
                 const std::string &output_path) const {
    if (!mesh_.containsPatch(patch_name)) {
      HERMES_LOG_ERROR("invalid patch id");
      return;
    }
    auto patch = mesh_.patch(patch_name);
    HERMES_LOG("Exporting patch {} to obj file {}", patch_name, output_path);

    // create map of vertices
    // global index -> obj index
    std::unordered_map<size_t, size_t> vertex_index_map;
    std::vector<hermes::point3> obj_vertices;
    // iterate over patch faces
    for (size_t i = 0; i < patch.size; ++i) {
      auto global_face_id = i + patch.start;
      // iterate face vertices
      for (auto v : mesh_.faces[global_face_id])
        if (!vertex_index_map.count(v)) {
          vertex_index_map[v] = obj_vertices.size();
          obj_vertices.emplace_back(mesh_.vertices[v]);
        }
    }

    // write obj
    hermes::Str content;
    for (const auto &vertex : obj_vertices)
      content.appendLine("v ", vertex.x, " ", vertex.y, " ", vertex.z);
    for (size_t i = 0; i < patch.size; ++i) {
      auto global_face_id = i + patch.start;
      content.append("f ");
      for (auto v : mesh_.faces[global_face_id])
        content.append(" ", vertex_index_map[v] + 1);
      content.appendLine();
    }
    hermes::FileSystem::writeFile(output_path, content.str());
  }
  ///
  /// \return Number of poly mesh patches
  [[nodiscard]] size_t patchCount() const { return mesh_.patches.size(); }
  /// Retrieves the id for a given patch name
  /// \param patch_name
  /// \return -1 if patch not found
  [[nodiscard]] int patchId(const std::string &patch_name) const {
    for (int i = 0; i < mesh_.patches.size(); ++i)
      if (patch_name == mesh_.patches[i].name)
        return i;
    return -1;
  }
  /// Patch information containing (name, type, group, start, size)
  /// \param patch_id
  /// \return
  [[nodiscard]] PolyMesh::Patch patchData(int patch_id) const {
    return mesh_.patches[patch_id];
  }
  /// Computes the face's center in a patch for a given face in local patch
  /// index \param patch_id \param local_face_id \return center coordinates
  /// (x,y,z)
  [[nodiscard]] py::array_t<f32> patchFaceCenter(size_t patch_id,
                                                 size_t local_face_id) const {
    auto result = py::array_t<f32>(3);
    auto r = result.mutable_unchecked<1>();
    auto c = mesh_.faceCenter(mesh_.patches[patch_id].start + local_face_id);
    for (int i = 0; i < 3; ++i)
      r[i] = c[i];
    return result;
  }
  /// Intersects a poly mesh patch with a given plane
  /// \param patch_name
  /// \param normal plane normal
  /// \param offset  plane offset
  /// \param sort_direction direction in which the intersecting faces must be
  /// sorted \return List of intersecting patch faces
  py::array_t<size_t> intersectPatch(const std::string &patch_name,
                                     const py::array_t<float> &normal,
                                     const py::array_t<float> &offset,
                                     const py::array_t<float> &sort_direction) {
    auto patch_id = patchId(patch_name);
    if (patch_id < 0)
      return {};
    auto n = normal.unchecked<1>();
    auto o = offset.unchecked<1>();
    auto s = sort_direction.unchecked<1>();
    auto sd = hermes::vec3(s[0], s[1], s[2]);
    hermes::Plane plane({n[0], n[1], n[2]}, {o[0], o[1], o[2]});
    std::vector<size_t> face_ids;
    // iterate through patch faces and collect intersected ones
    mesh_.iteratePatchFaces(patch_name,
                            [&](size_t face_id, size_t local_face_id) {
                              if (mesh_.faceIntersectsPlane(plane, face_id))
                                face_ids.emplace_back(face_id);
                            });
    std::sort(face_ids.begin(), face_ids.end(), [&](size_t a, size_t b) {
      return hermes::dot((hermes::vec3)mesh_.faceCenter(a), sd) <
             hermes::dot((hermes::vec3)mesh_.faceCenter(b), sd);
    });
    auto result = py::array_t<size_t>(face_ids.size());
    auto r = result.mutable_unchecked<1>();
    for (size_t i = 0; i < face_ids.size(); ++i)
      r[i] = face_ids[i] - mesh_.patches[patch_id].start;
    return result;
  }
  // *******************************************************************************************************************
  //                                                                                                            DEBUG
  // *******************************************************************************************************************
  /// Retrieves simulation debugging data
  /// \return
  [[nodiscard]] std::map<std::string, std::vector<double>> debugData() const {
    auto file_path = sim_path_.root() / ("debug_data/variables_" +
                                         sim_path_.currentFramePath().name());
    if (!file_path.exists())
      return {};
    auto lines = hermes::FileSystem::readLines(file_path);
    if (lines.empty())
      return {};
    std::map<std::string, std::vector<double>> debug_data;
    for (auto &line : lines) {
      auto vs = hermes::Str::split(hermes::Str::strip(line, "\n"));
      HERMES_CHECK_EXP(vs.size() == 21);
      debug_data["global_face_id"].emplace_back(std::stoi(vs[0]));
      debug_data["local_face_id"].emplace_back(std::stoi(vs[1]));
      debug_data["fluid_velocity.x"].emplace_back(std::stof(vs[2]));
      debug_data["fluid_velocity.y"].emplace_back(std::stof(vs[3]));
      debug_data["fluid_velocity.z"].emplace_back(std::stof(vs[4]));
      debug_data["dsl_velocity.x"].emplace_back(std::stof(vs[5]));
      debug_data["dsl_velocity.y"].emplace_back(std::stof(vs[6]));
      debug_data["dsl_velocity.z"].emplace_back(std::stof(vs[7]));
      debug_data["dsl_area"].emplace_back(std::stod(vs[8]));
      debug_data["snow_cover_mass"].emplace_back(std::stod(vs[9]));
      debug_data["dm"].emplace_back(std::stod(vs[10]));
      debug_data["entrainment_mass"].emplace_back(std::stod(vs[11]));
      debug_data["front_weight"].emplace_back(std::stod(vs[12]));
      debug_data["dsl_velocity_mag"].emplace_back(std::stod(vs[13]));
      debug_data["dsl_height"].emplace_back(std::stod(vs[14]));
      debug_data["ce"].emplace_back(std::stod(vs[15]));
      debug_data["we"].emplace_back(std::stod(vs[16]));
      debug_data["slope_angle"].emplace_back(std::stod(vs[17]));
      debug_data["entrained_depth"].emplace_back(std::stod(vs[18]));
      debug_data["acc_entrained_depth"].emplace_back(std::stod(vs[19]));
      debug_data["snow_cover_height"].emplace_back(std::stod(vs[20]));
    }
    return debug_data;
  }

private:
  // sim
  SimPath sim_path_;
  // mesh
  PolyMesh mesh_;
  // data fields
  std::vector<std::string> scalar_fields_;
  std::vector<std::string> vector_fields_;
};

// *********************************************************************************************************************
//                                                                                                       MeshUtils_py
// *********************************************************************************************************************
/// \brief OpenFOAM's PolyMesh wrapper
class PolyMesh_py {
public:
  // *******************************************************************************************************************
  //                                                                                                   STATIC METHODS
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  //                                                                                                 FRIEND FUNCTIONS
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  PolyMesh_py(bool verbose = false) { initHermes(verbose); }
  //                                                                                                       assignment
  // *******************************************************************************************************************
  //                                                                                                        OPERATORS
  // *******************************************************************************************************************
  //                                                                                                       assignment
  //                                                                                                       arithmetic
  //                                                                                                          boolean
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  void load(const std::string &msh_file_path) {
    psa_anim::MshReader reader;
    auto msh = reader.loadFromFile(msh_file_path);
    if (!msh) {
      HERMES_LOG_ERROR("failed to load the msh file {}", msh_file_path);
      return;
    }
    std::vector<hermes::point3> vertices;
    std::vector<std::vector<size_t>> faces;
    std::vector<std::vector<size_t>> cells;
    msh->retrieveMesh(vertices, faces, cells);
    //    mesh_ = PolyMesh::from(vertices, faces, cells);
  }
  // *******************************************************************************************************************
  //                                                                                                    PUBLIC FIELDS
  // *******************************************************************************************************************
private:
  PolyMesh mesh_;
};

PYBIND11_MAKE_OPAQUE(std::vector<std::string>);

PYBIND11_MODULE(psa_anim_py, m) {
  py::bind_vector<std::vector<std::string>>(m, "VectorString");
  m.doc() = "OpenFoam PolyMesh & Sim files wrapper";

  m.def("convert_to_quad_mesh", &stl2quad);

  m.def("terrain_base", &terrainBase);

  py::class_<BlockMeshTerrain_py>(m, "BlockMeshDesc")
      .def(py::init(
          [](bool verbose) { return new BlockMeshTerrain_py(verbose); }))
      .def("addPatchDirection", &BlockMeshTerrain_py::addPatchDirection,
           "define patch normal direction")
      .def("setBoundaryType", &BlockMeshTerrain_py::setBoundaryType,
           "define patch type")
      .def("setSlopeDirection", &BlockMeshTerrain_py::setSlopeDirection,
           "define XY slope direction")
      .def("setSlopePoint", &BlockMeshTerrain_py::setSlopePoint,
           "define XY slope point")
      .def("setHeight", &BlockMeshTerrain_py::setHeight, "define Z height")
      .def("setTopIsInclined", &BlockMeshTerrain_py::setTopIsInclined,
           "choose if top patch is inclined")
      .def("setFaceMatching", &BlockMeshTerrain_py::setFaceMatching,
           "choose if block faces should match")
      .def("setBlockResolution", &BlockMeshTerrain_py::setBlockResolution,
           "choose block resolution")
      .def("setBlockGrading", &BlockMeshTerrain_py::setGrading,
           "choose block grading")
      .def("loadOBJ", &BlockMeshTerrain_py::loadOBJ,
           "load terrain surface (quad mesh)")
      .def("save", &BlockMeshTerrain_py::save, "exports BlockMeshDict");

  py::class_<PolyMesh_py>(m, "OFPolyMesh")
      .def(py::init([](bool verbose) { return new PolyMesh_py(verbose); }))
      .def("loadFile", &PolyMesh_py::load, "load poly mesh from file (.msh)");

  py::class_<PSLDebugData>(m, "PSLDebugData")
      .def_readonly("global_face_id", &PSLDebugData::global_face_id)
      .def_readonly("local_face_id", &PSLDebugData::local_face_id)
      .def_readonly("dsl_area", &PSLDebugData::dsl_area)
      .def_readonly("snow_cover_mass", &PSLDebugData::snow_cover_mass)
      .def_readonly("dm", &PSLDebugData::dm)
      .def_readonly("entrainment_mass", &PSLDebugData::entrainment_mass)
      .def_readonly("front_weight", &PSLDebugData::front_weight)
      .def_readonly("dsl_velocity_mag", &PSLDebugData::dsl_velocity_mag)
      .def_readonly("dsl_height", &PSLDebugData::dsl_height)
      .def_readonly("we", &PSLDebugData::we)
      .def_readonly("slope_angle", &PSLDebugData::slope_angle);

  py::class_<DSL_py>(m, "DSL")
      .def(py::init([](bool verbose) { return new DSL_py(verbose); }))
      .def("patchName", &DSL_py::patchName, "patch name")
      .def("size", &DSL_py::size, "face count")
      .def("preparePSLInput", &DSL_py::preparePSLInput, "write PSL input")
      .def("extractSurface", &DSL_py::extractSurface,
           "extract polygonal surface mesh");

  py::class_<ProceduralDSL_py>(m, "ProceduralDSL")
      .def(py::init([](bool verbose) { return new ProceduralDSL_py(verbose); }))
      .def("patchName", &ProceduralDSL_py::patchName, "patch name")
      .def("size", &ProceduralDSL_py::size, "face count")
      .def("setDuration", &ProceduralDSL_py::setDuration, "total duration time")
      .def("setWriteInterval", &ProceduralDSL_py::setWriteInterval,
           "write interval time")
      .def("addRotatingDisk", &ProceduralDSL_py::addRotatingDisk,
           "add rotating disk")
      .def("addMovingBar", &ProceduralDSL_py::addMovingBar, "add moving bar")
      .def("addSteadyFront", &ProceduralDSL_py::addSteadyFront,
           "add steady front")
      .def("write", &ProceduralDSL_py::write, "write data");

  py::class_<PolyMesh::Patch>(m, "OFPatchData")
      .def_readonly("name", &PolyMesh::Patch::name)
      .def_readonly("type", &PolyMesh::Patch::type)
      .def_readonly("group", &PolyMesh::Patch::group)
      .def_readonly("start", &PolyMesh::Patch::start)
      .def_readonly("size", &PolyMesh::Patch::size);

  py::class_<OF_py>(m, "OFSim")
      .def(py::init([](bool verbose) { return new OF_py(verbose); }))
      .def("setSimPath", &OF_py::setSimPath, "open sim folder")
      .def("setFramesPath", &OF_py::setFramesPath, "open sim franes folder")
      .def("setCurrentFrame", &OF_py::setCurrentFrame, "set frame")
      .def("currentFramePath", &OF_py::currentFramePath,
           "get current frame path")
      .def("scalarFields", &OF_py::scalarFields, "list of scalar fields")
      .def("scalarFieldData", &OF_py::scalarFieldData,
           "list of scalar field values for a given patch")
      .def("debugData", &OF_py::debugData, "list of scalar debug values")
      .def("patchId", &OF_py::patchId, "find patch id from name")
      .def("patchCount", &OF_py::patchCount, "patch count")
      .def("patch", &OF_py::patchData, "patch info")
      .def("DSL", &OF_py::DSL, "DSL patch")
      .def("proceduralDSL", &OF_py::proceduralDSL, "procedural DSL patch")
      .def("patchFaceCenter", &OF_py::patchFaceCenter,
           "patch face center coordinates")
      .def("intersectPatch", &OF_py::intersectPatch,
           "find patch faces intersecting plane")
      .def("frameCount", &OF_py::frameCount, "simulation frame count")
      .def("exportOBJ", &OF_py::exportOBJ,
           "exports patch geometry into obj format");
}
