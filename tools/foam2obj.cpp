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
///\file foam2vdb.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2023-12-11
///
///\brief

#include <hermes/common/arg_parser.h>
#include <hermes/common/file_system.h>
#include <hermes/geometry/bbox.h>
#include <hermes/logging/logging.h>
#include <mesh_utils.h>
#include <openfoam_poly_mesh.h>

void writeObj(const PolyMesh &mesh, const std::set<size_t> &faces,
              bool wireframe, const hermes::Path &file) {
  hermes::Str content;
  std::unordered_map<size_t, size_t> vertex_map;
  for (auto f : faces) {
    // write vertices
    for (auto v : mesh.faces[f]) {
      if (!vertex_map.count(v)) {
        vertex_map[v] = vertex_map.size() + 1;
        content = content << "v " << mesh.vertices[v].x << " ";
        content = content << mesh.vertices[v].y << " " << mesh.vertices[v].z
                          << "\n";
      }
    }
    // create poly
    content = content << (wireframe ? "l " : "f ");
    for (auto v : mesh.faces[f])
      content = content << vertex_map[v] << " ";
    content.appendLine();
  }
  HERMES_LOG_VARIABLE(file.fullName());
  hermes::FileSystem::writeFile(file, content.str());
}

void writeBoundaryObj(const PolyMesh &mesh, const std::set<size_t> &region,
                      const hermes::Path &file) {
  auto topology = psa_anim::computeFaceTopology(mesh.faces, mesh.cells, region);
  std::set<size_t> faces;
  for (auto face : topology) {
    if (face.second.second < 0) {
      faces.insert(face.first);
    }
  }
  writeObj(mesh, faces, false, file);
}

int main(int argc, const char **argv) {
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  hermes::ArgParser parser("foam2obj", "converts a OpenFOAM mesh to obj file.");
  parser.addArgument("-i", "OpenFOAM simulation path", true);
  parser.addArgument("-o", "OpenVDB output path", true);
  parser.addArgument(
      "--bbox", "Convert only the cells that intersect with the bounding box");
  parser.addArgument("--boundary", "Output region boundary faces");
  parser.parse(argc, argv, true);

  hermes::Path path(parser.get<std::string>("-i"));
  hermes::Path output_path(parser.get<std::string>("-o"));
  bool boundary = parser.check("--boundary");

  hermes::bbox3 bbox;
  if (parser.check("--bbox")) {
    std::vector<float> values = parser.getList<float>("--bbox");
    bbox = hermes::make_union(bbox,
                              hermes::point3(values[0], values[1], values[2]));
    bbox = hermes::make_union(bbox,
                              hermes::point3(values[3], values[4], values[5]));
  }

  PolyMesh mesh;
  mesh.load(path / "constant/polyMesh");

  mesh.computeCellBounds();
  size_t cell_count = mesh.cells.size();

  HERMES_LOG_VARIABLE(mesh.faces.size());

  std::set<size_t> output_faces;
  std::set<size_t> output_cells;
  for (size_t cell_id = 0; cell_id < cell_count; ++cell_id) {
    if (parser.check("--bbox") &&
        !hermes::overlaps(bbox, mesh.cell_bounds[cell_id]))
      continue;
    output_cells.insert(cell_id);
    for (auto f : mesh.cells[cell_id])
      output_faces.insert(f);
  }

  HERMES_LOG_VARIABLE(output_faces.size());

  if (boundary)
    writeBoundaryObj(mesh, output_cells, output_path);
  else
    writeObj(mesh, output_faces, true, output_path);

  return 0;
}
