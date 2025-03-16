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
///\file poly_mesh.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-03-11
///
///\brief

#ifndef PSA_ANIM_TOOLS_FOAM2VDB_POLY_MESH_H
#define PSA_ANIM_TOOLS_FOAM2VDB_POLY_MESH_H

#include "hermes/geometry/queries.h"
#include "hermes/logging/logging.h"
#include <hermes/geometry/bbox.h>
#include <hermes/geometry/plane.h>
#include <hermes/geometry/point.h>
#include <mesh_utils.h>
#include <openfoam_parser.h>
#include <set>
#include <vector>

// *********************************************************************************************************************
//                                                                                                           PolyMesh
// *********************************************************************************************************************
/// \brief Holds an OpenFOAM's PolyMesh data structure.
/// Based on https://doc.cfd.direct/openfoam/user-guide-v10/mesh-description ,
/// the PolyMesh contains:
///
/// Vertices:
///     (x, y, z) vertex positions
///
/// Faces:
///  A face is an ordered list of points, where a point is referred to by its
///  label. The ordering of point labels in a face is such that each two
///  neighbouring points are connected by an edge.
///
///  Internal faces
///     Those faces that connect two cells (and it can never be more than two).
///     For each internal face, the ordering of the point labels is such that
///     the face normal points into the cell with the larger label, i.e. for
///     cells 2 and 5, the normal points into 5;
///  Boundary faces
///     Those belonging to one cell since they coincide with the boundary of the
///     domain. A boundary face is therefore addressed by one cell(only) and a
///     boundary patch. The ordering of the point labels is such that the face
///     normal points outside of the computational domain.
///
/// Cells:
///  A cell is a list of faces in arbitrary order. Cells must have the
///  properties listed below. Contiguous
///     The cells must completely cover the computational domain and must not
///     overlap one another.
///  Convex
///     Every cell must be convex and its cell centre inside the cell.
///  Closed
///     Every cell must be closed, both geometrically and topologically where:
///         - geometrical closedness requires that when all face area vectors
///         are oriented to point outwards of the
///           cell, their sum should equal the zero vector to machine accuracy;
///         - topological closedness requires that all the edges in a cell are
///         used by exactly two faces of the cell
///           in question.
///
/// Boundary:
///  A boundary is a list of patches, each of which is associated with a
///  boundary condition. A patch is a list of face labels which clearly must
///  contain only boundary faces and no internal faces. The boundary is required
///  to be closed, i.e. the sum all boundary face area vectors equates to zero
///  to machine tolerance.
/// -----------------------------------------------------------------------------------------------------------------
/// The PolyMesh is described by the following data arrays:
///  [points]
///     a list of vectors describing the cell vertices, where the first vector
///     in the list represents vertex 0, the second vector represents vertex 1,
///     etc.;
///  [faces]
///     a list of faces, each face being a list of indices to vertices in the
///     points list, where again, the first entry in the list represents face 0,
///     etc.;
///  [owner]
///     a list of owner cell labels, the index of entry relating directly to the
///     index of the face, so that the first entry in the list is the owner
///     label for face 0, the second entry is the owner label for face 1, etc;
///  [neighbour]
///     a list of neighbour cell labels (neighbours are defined only for
///     internal faces);
///  [boundary]
///     a list of patches, containing a dictionary entry for each patch,
///     declared using the patch name
struct PolyMesh {
  // *******************************************************************************************************************
  //                                                                                                  PolyMesh::Patch
  // *******************************************************************************************************************
  struct Patch {
    std::string name;  //!< patch name identifier
    std::string type;  //!< patch type: wall, empty, patch, ...
    std::string group; //!< patch group id
    size_t start{};    //!< global face index start
    size_t size{};     //!< patch face count
  };
  // *******************************************************************************************************************
  //                                                                                                   STATIC METHODS
  // *******************************************************************************************************************
  /// \brief Constructs a poly mesh from a list of vertices, faces and cells
  /// \note This function assumes:
  /// \note - Cells are closed and convex.
  /// \note - Face's vertex indices are ordered (the normal direction is defined
  /// by the right-hand rule). \note - Faces are sorted by boundary patches.
  /// \note - Boundary face's normals point towards the outside.
  /// \note - Boundary faces come sorted at the end.
  /// \param _vertices
  /// \param _faces
  /// \param _cells
  /// \return
  static hermes::Result<PolyMesh>
  from(const std::vector<hermes::point3> &_vertices,
       const std::vector<std::vector<size_t>> &_faces,
       const std::vector<std::vector<size_t>> &_cells,
       const std::vector<Patch> &_patches) {
    // return object
    PolyMesh mesh;

    // copy data
    mesh.patches = _patches;
    mesh.vertices = _vertices;
    mesh.cells = _cells;
    mesh.faces = _faces;
    mesh.cell_centers = psa_anim::cellCenters(_vertices, _faces, _cells);

    // The first restriction is about face ordering:
    auto face_cell_map = psa_anim::computeFaceTopology(mesh.faces, mesh.cells);
    size_t internal_face_count = 0;
    for (auto &face_cells : face_cell_map) {
      auto first_cell_id = face_cells.second.first;
      auto second_cell_id = face_cells.second.second;
      auto face_id = face_cells.first;
      auto cell_center = mesh.cell_centers[face_cells.second.first];
      auto normal =
          psa_anim::faceNormal(mesh.vertices, mesh.faces[face_cells.first]);
      auto face_center =
          psa_anim::faceCenter(mesh.vertices, mesh.faces[face_id]);
      auto center_direction = cell_center - face_center;
      bool normal_points_first_cell_center =
          hermes::dot(center_direction, normal) > 0;
      bool flip_face = false;
      if (second_cell_id < 0) {
        // for boundary faces, the normal should point towards the outside of
        // the domain sanity check as input is assumed to fulfill this
        // requirement face should point away from cell center
        flip_face = normal_points_first_cell_center;
      } else {
        internal_face_count++;
        HERMES_ASSERT(first_cell_id > second_cell_id);
        // for internal faces, the normal should point towards the cell with
        // larger index and we expect the first index to be the larger
        flip_face = !normal_points_first_cell_center;
      }
      if (flip_face) {
        auto face_copy = mesh.faces[face_id];
        for (size_t i = 0; i < face_copy.size(); ++i)
          face_copy[i] = mesh.faces[face_id][face_copy.size() - i - 1];
      }
    }

    // cell faces should be ordered by the cell index (upper triangular order)
    auto otherCell = [&](size_t cell_id, size_t face_id) -> i64 {
      auto it = face_cell_map.find(face_id);
      HERMES_ASSERT(it != face_cell_map.end());
      HERMES_ASSERT(it->second.first == cell_id ||
                    it->second.second == cell_id);
      return it->second.first == cell_id ? it->second.second : it->second.first;
    };
    for (size_t cell_id = 0; cell_id < mesh.cells.size(); ++cell_id) {
      std::sort(mesh.cells[cell_id].begin(), mesh.cells[cell_id].end(),
                [&](size_t face_a, size_t face_b) {
                  return otherCell(cell_id, face_a) <
                         otherCell(cell_id, face_b);
                });
    }

    // owner cell labels are the first index in the face_cell_map
    mesh.owner.resize(mesh.faces.size());
    // neighbour cell labels are the second index in the face_cell_map
    mesh.neighbour.resize(internal_face_count);

    for (auto &face_cells : face_cell_map) {
      mesh.owner[face_cells.first] = face_cells.second.first;
      if (face_cells.first < internal_face_count)
        mesh.neighbour[face_cells.first] = face_cells.second.second;
    }

    return hermes::Result<PolyMesh>(std::move(mesh));
  }
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  //                                                                                                               IO
  /// Loads points from openfoam dict file.
  /// \param path points file path
  void readPoints(const hermes::Path &path) {
    OpenFoamDict dict(path);
    if (dict.nodes().count("list"))
      vertices =
          OpenFoamDict::parseValuesFrom<hermes::point3>(dict["list"].value);
    HERMES_LOG("Read {} points.", vertices.size());
  }
  /// Loads faces from openfoam dict file.
  /// \param path faces file path
  void readFaces(const hermes::Path &path) {
    OpenFoamDict dict(path);
    if (dict.nodes().count("list"))
      OpenFoamDict::parseValueListsFrom<size_t>(dict["list"].value, faces);
    HERMES_LOG("Read {} faces.", faces.size());
  }
  /// Loads boundary from openfoam dict file.
  /// \param path boundary file path
  void readBoundary(const hermes::Path &path) {
    OpenFoamDict dict(path);
    if (dict.nodes().count("list")) {
      auto nodes = OpenFoamDict::parseValuesFrom<OpenFoamDict::DictNode>(
          dict["list"].value);
      for (const auto &node : nodes) {
        Patch patch;
        patch.name = node.name;
        patch.type = node.fields.count("type") ? node["type"].value : "";
        patch.group =
            node.fields.count("inGroups") ? node["inGroups"].value : "";
        patch.size =
            node.fields.count("nFaces") ? std::stoul(node["nFaces"].value) : 0;
        patch.start = node.fields.count("startFace")
                          ? std::stoul(node["startFace"].value)
                          : 0;
        patches.emplace_back(patch);
      }
    }
    HERMES_LOG("Read {} patches.", patches.size());
  }
  /// Loads owners from openfoam dict file.
  /// \param path owner file path
  void readOwners(const hermes::Path &path) {
    OpenFoamDict dict(path);
    if (dict.nodes().count("list"))
      owner = OpenFoamDict::parseValuesFrom<size_t>(dict["list"].value);
    HERMES_LOG("Read {} owners.", owner.size());
  }
  ///
  /// \param path
  void readNeighbours(const hermes::Path &path) {
    OpenFoamDict dict(path);
    if (dict.nodes().count("list"))
      neighbour = OpenFoamDict::parseValuesFrom<i64>(dict["list"].value);
    HERMES_LOG("Read {} neighbours.", neighbour.size());
  }
  /// Loads the entire poly mesh structure from a constant/polyMesh folder
  /// structure. \param path polyMesh folder path
  void load(const hermes::Path &path) {
    HERMES_LOG_VARIABLE(path);
    readPoints(path / "points");
    readFaces(path / "faces");
    readBoundary(path / "boundary");
    readOwners(path / "owner");
    readNeighbours(path / "neighbour");
    constructCells();
    computeCellCenters();
    computeCellBounds();
    HERMES_LOG_VARIABLE(faces.size());
    HERMES_LOG_VARIABLE(cells.size());
    HERMES_LOG_VARIABLE(vertices.size());
    HERMES_LOG_VARIABLE(owner.size());
    HERMES_LOG_VARIABLE(neighbour.size());
  }
  /// \brief Saves the poly mesh into the OpenFOAM's folder structure
  /// \note This function creates the following dict files: boundary, faces,
  /// neighbour, owner and points. \param path path to the polyMesh folder
  /// (destination of the files)
  void save(const hermes::Path &path) {
    if (!path.exists() || !path.isDirectory()) {
      HERMES_LOG_ERROR("invalid path.");
      return;
    }
    static const char *header_template =
        "/*--------------------------------*- C++ "
        "-*-----------------------------------*\\\n"
        "| =========                 |                                         "
        "          |\n"
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox "
        "          |\n"
        "|  \\\\    /   O peration     | Version:  2012                        "
        "          |\n"
        "|   \\\\  /    A nd           | Website:  www.openfoam.com            "
        "          |\n"
        "|    \\\\/     M anipulation  |                                       "
        "          |\n"
        "\\*-------------------------------------------------------------------"
        "---------*/\n"
        "FoamFile\n"
        "{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        "    class       <<<CLASS>>>;\n"
        "    <<<NOTE>>>\n"
        "    location    \"constant/polyMesh\";\n"
        "    object      <<<OBJECT>>>;\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
        "* * * * //\n";

    // compute boundary faces
    size_t boundary_face_count = 0;
    for (const auto &p : patches)
      boundary_face_count += p.size;
    auto internal_face_count = faces.size() - boundary_face_count;
    auto note = hermes::Str::format(
        "note     \"nPoints:{}  nCells:{}  nFaces:{}  nInternalFaces:{}\";\n",
        vertices.size(), cells.size(), faces.size(), internal_face_count);

    // store points
    {
      auto header = hermes::Str::regex::replace(header_template, "<<<CLASS>>>",
                                                "vectorField");
      header = hermes::Str::regex::replace(header, "<<<OBJECT>>>", "points");
      header = hermes::Str::regex::replace(header, "<<<NOTE>>>", "");
      hermes::Str content = header;
      content.appendLine(vertices.size());
      content.appendLine("(");
      for (const auto &point : vertices)
        content.appendLine("(", point.x, " ", point.y, " ", point.z, ")");
      content.appendLine(")");
      hermes::FileSystem::writeFile(path / "points", content.str());
    }

    // store faces
    {
      auto header = hermes::Str::regex::replace(header_template, "<<<CLASS>>>",
                                                "faceList");
      header = hermes::Str::regex::replace(header, "<<<OBJECT>>>", "faces");
      header = hermes::Str::regex::replace(header, "<<<NOTE>>>", "");
      hermes::Str content = header;
      content.appendLine(faces.size());
      content.appendLine("(");
      for (const auto &face : faces)
        content.appendLine(face.size(), "(", hermes::Str::join(face, " "), ")");
      content.appendLine(")");
      hermes::FileSystem::writeFile(path / "faces", content.str());
    }

    // store owner
    {
      auto header = hermes::Str::regex::replace(header_template, "<<<CLASS>>>",
                                                "labelList");
      header = hermes::Str::regex::replace(header, "<<<OBJECT>>>", "owner");
      header = hermes::Str::regex::replace(header, "<<<NOTE>>>", note);
      hermes::Str content = header;
      content.appendLine(faces.size());
      content.appendLine("(");
      for (size_t i = 0; i < faces.size(); ++i)
        content.appendLine(owner[i]);
      content.appendLine(")");
      hermes::FileSystem::writeFile(path / "owner", content.str());
    }

    // store neighbour
    {
      auto header = hermes::Str::regex::replace(header_template, "<<<CLASS>>>",
                                                "labelList");
      header = hermes::Str::regex::replace(header, "<<<OBJECT>>>", "neighbour");
      header = hermes::Str::regex::replace(header, "<<<NOTE>>>", note);
      hermes::Str content = header;
      content.appendLine(internal_face_count);
      content.appendLine("(");
      for (size_t i = 0; i < internal_face_count; ++i) {
        HERMES_ASSERT(neighbour[i] >= 0);
        content.appendLine(neighbour[i]);
      }
      content.appendLine(")");
      hermes::FileSystem::writeFile(path / "neighbour", content.str());
    }

    // store boundary
    {
      auto header = hermes::Str::regex::replace(header_template, "<<<CLASS>>>",
                                                "polyBoundaryMesh");
      header = hermes::Str::regex::replace(header, "<<<OBJECT>>>", "boundary");
      header = hermes::Str::regex::replace(header, "<<<NOTE>>>", "");
      hermes::Str content = header;
      content.appendLine(patches.size());
      content.appendLine("(");
      for (const auto &p : patches) {
        content.appendLine("\t", p.name);
        content.appendLine("\t{");
        content.appendLine("\t\ttype         ", p.type, ";");
        content.appendLine("\t\tnFaces       ", p.size, ";");
        content.appendLine("\t\tstartFace    ", p.start, ";");
        content.appendLine("\t}");
      }
      content.appendLine(")");
      hermes::FileSystem::writeFile(path / "boundary", content.str());
    }
  }
  //                                                                                                           access
  /// Retrieve patch by patch name
  /// \param name patch identifier
  /// \return patch object or empty object if not found
  [[nodiscard]] Patch patch(const std::string &name) const {
    for (const auto &p : patches)
      if (p.name == name)
        return p;
    HERMES_LOG_ERROR("patch not found!");
    return {};
  }
  /// Checks if patch with given name exists
  /// \param name patch identifier
  /// \return true if patch exists
  [[nodiscard]] bool containsPatch(const std::string &name) const {
    return std::any_of(patches.begin(), patches.end(),
                       [&](const Patch &p) { return p.name == name; });
    for (const auto &p : patches)
      if (p.name == name)
        return true;
    return false;
  }
  ///
  /// \param name
  /// \return
  [[nodiscard]] int patchId(const std::string &name) const {
    for (size_t i = 0; i < patches.size(); ++i)
      if (patches[i].name == name)
        return i;
    return -1;
  }
  //                                                                                                         geometry
  bool cellContainsPoint(size_t cell_id, const hermes::point3 &p) const {
    if (cell_id > cells.size()) {
      HERMES_LOG_ERROR("invalid cell id {}", cell_id);
      return false;
    }

    auto center = cell_centers[cell_id];
    for (auto f : cells[cell_id]) {
      auto face_normal = faceNormal(f);
      auto face_center = faceCenter(f);
      auto dc = hermes::dot(face_normal, center - face_center);
      if (dc < 0) {
        face_normal *= -1;
        dc *= -1;
      }
      auto dp = hermes::dot(face_normal, p - face_center);
      if (dc * dp < 0)
        return false;
    }
    return true;

    // shoot ray
    hermes::ray3 r(p, hermes::vec3(1, 0, 0));
    size_t intersection_count = 0;
    for (auto f : cells[cell_id]) {
      // face vertices
      auto vertices = faceVertices(f);
      if (psa_anim::polyFaceIntersect(vertices, r))
        intersection_count++;
    }
    return intersection_count % 2;
  }

  [[nodiscard]] size_t faceSize(size_t i) const {
    HERMES_ASSERT(i < faces.size());
    return faces[i].size();
  }

  [[nodiscard]] std::vector<hermes::point3> cellVertices(size_t i) const {
    std::vector<hermes::point3> cell_vertices;
    std::set<size_t> s;
    for (const auto &f : cells[i])
      for (const auto &v : faces[f])
        s.insert(v);
    for (auto v : s)
      cell_vertices.emplace_back(vertices[v]);
    return cell_vertices;
  }

  [[nodiscard]] hermes::vec3 faceNormal(size_t global_face_id) const {
    auto vs = faceVertices(global_face_id);
    auto a = vs[1] - vs[0];
    auto b = vs[2] - vs[0];
    return hermes::normalize(hermes::cross(a, b));
    return psa_anim::faceNormal(vertices, faces[global_face_id]);
  }

  //  [[nodiscard]] hermes::point3 cellCenter(size_t i) const {
  // iterate over cell faces
  //    auto cell_vertices = cellVertices(i);
  //    hermes::point3 cell_center;
  //    HERMES_CHECK_EXP(!cell_vertices.empty());
  //    for (auto v:cell_vertices)
  //      cell_center += (hermes::vec3) v;
  //    return cell_center / cell_vertices.size();
  //  }

  [[nodiscard]] hermes::bbox3 cellBounds(size_t i) const {
    // iterate over cell faces
    auto cell_vertices = cellVertices(i);
    hermes::bbox3 cell_bb;
    HERMES_CHECK_EXP(!cell_vertices.empty());
    for (auto v : cell_vertices)
      cell_bb = hermes::make_union(cell_bb, v);
    return cell_bb;
  }

  [[nodiscard]] hermes::point3 faceCenter(size_t i) const {
    if (i >= faces.size())
      HERMES_LOG_ERROR("face id out of bounds!");
    auto n = faces[i].size();
    hermes::point3 center;
    for (size_t j = 0; j < n; ++j)
      center += (hermes::vec3)vertices[faces[i][j]];
    center /= n;
    return center;
  }

  void iteratePatchFaces(const std::string &name,
                         const std::function<void(size_t, size_t)> &f) const {
    for (auto &p : patches)
      if (p.name == name) {
        for (size_t i = p.start; i < p.start + p.size; ++i)
          f(i, i - p.start);
        return;
      }
    HERMES_LOG_ERROR("patch not found!");
  }

  [[nodiscard]] bool faceIntersectsPlane(const hermes::Plane &plane,
                                         size_t global_face_id) const {
    auto vs = faceVertices(global_face_id);
    int side_count[2] = {0, 0};
    for (const auto &v : vs)
      if (plane.onNormalSide(v))
        side_count[0]++;
      else
        side_count[1]++;
    return side_count[0] && side_count[1];
  }

  void constructCells() {
    HERMES_LOG("Construct cells.");
    // get number of cells
    size_t n = 0;
    for (auto o : owner)
      n = std::max(o, n);
    HERMES_LOG("... found {} max index cell.", n + 1);
    for (auto o : neighbour) {
      if (o >= 0)
        n = std::max((size_t)o, n);
    }
    HERMES_LOG("... found {} max index cell.", n + 1);
    cells.resize(n + 1);
    for (size_t i = 0; i < owner.size(); ++i) {
      HERMES_ASSERT(owner[i] < cells.size());
      cells[owner[i]].emplace_back(i);
    }
    // the neighbour list exists only for internal faces
    for (size_t i = 0; i < neighbour.size(); ++i)
      if (neighbour[i] >= 0) {
        if (neighbour[i] >= cells.size()) {
          HERMES_LOG_VARIABLE(neighbour[i]);
          HERMES_LOG_VARIABLE(cells.size());
        }
        HERMES_ASSERT(neighbour[i] < cells.size());
        cells[neighbour[i]].emplace_back(i);
      }
  }

  void computeCellCenters() {
    HERMES_LOG("Compute cell centers.");
    cell_centers = psa_anim::cellCenters(vertices, faces, cells);
    //    cell_centers = psa_anim::cellCenters(vertices, faces, cells);
    //    cell_centers.resize(cells.size());
    //    for (size_t cell_id = 0; cell_id < cells.size(); ++cell_id)
    //      cell_centers[cell_id] = cellCenter(cell_id);
  }

  void computeCellBounds() {
    HERMES_LOG("Compute cell bounds.");
    cell_bounds.resize(cells.size());
    for (size_t cell_id = 0; cell_id < cells.size(); ++cell_id)
      cell_bounds[cell_id] = cellBounds(cell_id);
  }
  /// Get vertices from face
  /// \param global_face_id
  /// \return list of face vertices
  [[nodiscard]] std::vector<hermes::point3>
  faceVertices(size_t global_face_id) const {
    std::vector<hermes::point3> face_vertices;
    for (auto v : faces[global_face_id])
      face_vertices.emplace_back(vertices[v]);
    return face_vertices;
  }
  /// Get vertex indices from face
  /// \param global_face_id
  /// \return list of face vertices
  [[nodiscard]] std::vector<size_t>
  faceVertexIndices(size_t global_face_id) const {
    std::vector<size_t> face_vertices;
    for (auto v : faces[global_face_id])
      face_vertices.emplace_back(v);
    return face_vertices;
  }

  bool isFacingForwards(size_t face_id, const hermes::point3 &p) {
    HERMES_CHECK_EXP(faces[face_id].size() > 2);
    // lazily checks against first triangle
    auto a = vertices[faces[face_id][0]];
    auto b = vertices[faces[face_id][1]];
    auto c = vertices[faces[face_id][2]];
    auto n = hermes::cross(b - a, c - b);
    auto d = p - a;
    return hermes::dot(n, d) > 0;
  }

  bool checkFaceOrientation() {
    if (cell_centers.size() != cells.size())
      computeCellCenters();
    // for each cell, check if faces are oriented accordingly to owner/neighbour
    // relationship
    for (size_t cell_id = 0; cell_id < cells.size(); ++cell_id)
      for (size_t face_id : cells[cell_id]) {
        bool is_owner = owner[face_id] == cell_id;
        bool faces_forwards_center =
            isFacingForwards(face_id, cell_centers[cell_id]);
        if ((is_owner && faces_forwards_center) ||
            (!is_owner && !faces_forwards_center)) {
          HERMES_LOG_WARNING(
              "face is owner: {} face_forwards_center: {} face id: {}",
              (int)is_owner, (int)faces_forwards_center, face_id);
          return false;
        }
      }
    return true;
  }

  bool checkCentersInsideCells() {
    for (size_t cell_id = 0; cell_id < cells.size(); ++cell_id)
      if (!isInsideCell(cell_centers[cell_id], cell_id))
        return false;
    return true;
  }

  bool isInsideCell(const hermes::point3 &p, size_t cell_id) {
    for (size_t face_id : cells[cell_id]) {
      bool is_owner = owner[face_id] == cell_id;
      bool faces_forwards_center = isFacingForwards(face_id, p);
      if ((is_owner && faces_forwards_center) ||
          (!is_owner && !faces_forwards_center))
        return false;
    }
    return true;
  }

  //                                                                                                            check
  ///
  /// \return
  [[nodiscard]] bool check() const {
    // mesh checks:
    // 0 - data checks
    // 1 - closed cells
    // 2 - face centers
    // 3 - face usage
    // 4 - owner / neighbour -> normal consistency
    // 5 - patches
    // 6 - cell sizes
    // 7 - cell centers
    // 8 - boundary faces at the end

    // 0 - sanity checks
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(!cells.empty(), false,
                                       "empty cell data");
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(!faces.empty(), false,
                                       "empty face data");
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(!vertices.empty(), false,
                                       "empty point data");
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(!owner.empty(), false,
                                       "empty owner data");
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(!neighbour.empty(), false,
                                       "empty neighbour data");
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(cell_centers.size() == cells.size(),
                                       false, "#cells != #cell_centers");

    // auxiliary data
    auto face_cells = psa_anim::computeFaceTopology(faces, cells);

    // 6 - cell sizes
    // from face topology, we should retrieve all cells
    std::unordered_map<size_t, size_t> cell_sizes;
    for (const auto &face_topo : face_cells) {
      cell_sizes[face_topo.second.first]++;
      cell_sizes[face_topo.second.second]++;
    }
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(cell_sizes.size() == cells.size() + 1,
                                       false, "empty cells");
    std::unordered_map<size_t, size_t> size_frequency;
    for (size_t cell_id = 0; cell_id < cells.size(); ++cell_id) {
      HERMES_LOG_AND_RETURN_VALUE_IF_NOT(cell_sizes.count(cell_id), false,
                                         "empty cells");
      size_frequency[cell_sizes[cell_id]]++;
    }
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(cell_sizes.count(-1), false,
                                       "no boundary cells");
    for (const auto &sf : size_frequency)
      HERMES_LOG("size {} -> frequency {}", sf.first, sf.second);

    // 1 - closed cells (all face edges must be shared by exactly two faces of
    // the cell)
    auto edgeKey = [](size_t a, size_t b) -> std::pair<size_t, size_t> {
      HERMES_ASSERT(a != b);
      return {std::min(a, b), std::max(a, b)};
    };
    for (const auto &cell_faces : cells) {
      // create edge usage map for this cell
      //                edge key, face count
      std::map<std::pair<size_t, size_t>, size_t> edge_usage;
      // register edges
      for (const auto &face_id : cell_faces) {
        HERMES_ASSERT(face_id < faces.size());
        size_t face_size = faces[face_id].size();
        for (size_t i = 0; i < face_size; i++)
          edge_usage[edgeKey(faces[face_id][i],
                             faces[face_id][(i + 1) % face_size])]++;
      }
      // all edges should be used exactly twice
      for (const auto &usage : edge_usage) {
        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(usage.second == 2, false,
                                           "edge usage in cell != 2");
      }
    }

    // 2 - face centers (face centers should be "inside" faces)
    for (const auto &face_vertices : faces) {
      auto center = psa_anim::faceCenter(vertices, face_vertices);
      HERMES_LOG_AND_RETURN_VALUE_IF_NOT(face_vertices.size() > 2, false,
                                         "#face_vertices < 3");
      bool common_side = false;
      for (size_t i = 0; i < face_vertices.size(); ++i) {
        auto a = i;
        auto b = (i + 1) % face_vertices.size();
        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(face_vertices[a] < vertices.size(),
                                           false, "invalid face vertex index");
        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(face_vertices[b] < vertices.size(),
                                           false, "invalid face vertex index");
        auto va = vertices[b] - vertices[a];
        auto vb = center - vertices[a];
        auto side = hermes::dot(va, vb) < 0;
        if (i == 0)
          common_side = side;
        //        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(common_side == side, false,
        //        "face center outside face");
      }
    }

    // 3 - face usage (all faces should be used by at most 2 cells, and no less
    // than 1 cell)
    std::unordered_map<size_t, size_t> face_usage;
    for (const auto &cell_faces : cells)
      for (auto face_id : cell_faces)
        face_usage[face_id]++;
    for (size_t i = 0; i < faces.size(); ++i) {
      HERMES_LOG_AND_RETURN_VALUE_IF_NOT(
          face_usage.count(i), false,
          "found an alone face (no neighbour cells)");
      HERMES_LOG_AND_RETURN_VALUE_IF_NOT(face_usage[i] > 0 && face_usage[i] < 3,
                                         false, "bad face sharing");
    }

    // 4 - owner / neighbour -> normal consistency
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(face_cells.size() == faces.size(), false,
                                       "bad topology");
    // check for topology
    for (const auto &face_t : face_cells)
      HERMES_LOG_AND_RETURN_VALUE_IF_NOT(
          face_t.second.first > face_t.second.second, false, "bad topology");
    // here we also count the number of boundary faces to match with the patch
    // data later
    size_t boundary_face_count = 0;
    size_t internal_face_count = 0;
    for (const auto &face_topology : face_cells) {
      auto face_id = face_topology.first;
      auto face_normal = psa_anim::faceNormal(vertices, faces[face_id]);
      auto face_center = psa_anim::faceCenter(vertices, faces[face_id]);
      auto cell_a = face_topology.second.first;
      auto cell_b = face_topology.second.second;
      HERMES_LOG_AND_RETURN_VALUE_IF_NOT(cell_a != cell_b, false,
                                         "bad face topology");
      HERMES_LOG_AND_RETURN_VALUE_IF_NOT(cell_a >= 0, false,
                                         "found face with no owner cell");
      HERMES_LOG_AND_RETURN_VALUE_IF_NOT(cell_a < cells.size(), false,
                                         "invalid face owner");
      auto center_a = cell_centers[cell_a];
      auto center_a_direction = center_a - face_center;
      auto points_a = hermes::dot(center_a_direction, face_normal) > 0;
      if (cell_b >= 0) {
        internal_face_count++;
        // Those faces that connect two cells (and it can never be more than
        // two). For each internal face, the ordering of the point labels is
        // such that the face normal points into the cell with the larger label.
        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(
            points_a, false, "inconsistent internal face normal direction");
      } else {
        boundary_face_count++;
        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(!points_a, false,
                                           "boundary normal inwards");
      }
    }

    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(owner.size() == faces.size(), false,
                                       "bad #owner");
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(neighbour.size() == internal_face_count,
                                       false, "bad #neighbour");

    for (const auto &face_topology : face_cells) {
      auto face_id = face_topology.first;
      auto cell_a = face_topology.second.first;
      auto cell_b = face_topology.second.second;
      //      if (face_topology.first < internal_face_count)
      //        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(cell_b == neighbour[face_id],
      //        false, "bad neighbour index");
      // check owner
      //      HERMES_LOG_AND_RETURN_VALUE_IF_NOT(cell_a == owner[face_id],
      //      false, "bad owner index");
    }

    // 5 - patches (patch sizes == boundary faces)
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(
        boundary_face_count + internal_face_count == faces.size(), false,
        "bad boundary/internal face count");
    auto patch_sum = 0;
    for (const auto &p : patches)
      patch_sum += p.size;
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(patch_sum == boundary_face_count, false,
                                       "bad patch face coverage");
    // check for overlapping patches
    for (size_t i = 0; i < patches.size(); ++i)
      for (size_t j = i + 1; j < patches.size(); ++j) {
        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(patches[i].start != patches[j].start,
                                           false, "overlapping patches");
        if (patches[i].start > patches[j].start) {
          HERMES_LOG_AND_RETURN_VALUE_IF_NOT(
              patches[i].start >= patches[j].start + patches[j].size, false,
              "overlapping patches");
        } else {
          HERMES_LOG_AND_RETURN_VALUE_IF_NOT(
              patches[j].start >= patches[i].start + patches[i].size, false,
              "overlapping patches");
        }
      }

    // 7 - cell centers (should lie inside the cell)
    for (size_t cell_id = 0; cell_id < cell_centers.size(); ++cell_id) {
      // for each face, check the normal
    }

    // 8 - boundary faces at the end
    // patch faces should be boundary and at the end
    for (const auto &p : patches)
      for (size_t i = p.start; i < p.start + p.size; ++i) {
        HERMES_ASSERT(face_cells.count(i));
        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(i >= internal_face_count, false,
                                           "patch face not at the end");
        HERMES_LOG_AND_RETURN_VALUE_IF_NOT(face_cells[i].second < 0, false,
                                           "patch face not a boundary face");
      }

    return true;
  }

  psa_anim::FaceMesh intersect(const hermes::Plane &plane) {

    psa_anim::FaceMesh mesh;
    // all intersection points
    // a -> b -> vertex index
    std::unordered_map<size_t, std::unordered_map<size_t, size_t>> edge_points;

    auto contains = [&](size_t a, size_t b) -> bool {
      auto it = edge_points.find(a);
      if (it != edge_points.end())
        return it->second.count(b);
      return false;
    };

    auto intersectEdge = [&](const hermes::point3 &a, const hermes::point3 &b,
                             hermes::point3 &c) -> bool {
      hermes::Ray3 ray(a, b - a);
      if (auto r = hermes::GeometricPredicates::intersect(plane, ray)) {
        if (plane.onNormalSide(a) == plane.onNormalSide(b))
          return false;
        c = ray(*r);
        return true;
      }
      return false;
    };

    for (const auto &cell : cells) {
      std::vector<size_t> face_vertices;
      for (auto face_id : cell) {
        std::vector<size_t> face_edge_vertices;
        for (size_t i = 0; i < faces[face_id].size(); ++i) {
          auto edge_a = faces[face_id][i];
          auto edge_b = faces[face_id][(i + 1) % faces[face_id].size()];
          size_t a = std::min(edge_a, edge_b);
          size_t b = std::max(edge_a, edge_b);
          if (contains(a, b)) {
            face_edge_vertices.emplace_back(edge_points[a][b]);
            continue;
          }
          hermes::point3 p;
          if (intersectEdge(vertices[a], vertices[b], p)) {
            edge_points[a][b] = mesh.vertices.size();
            face_edge_vertices.emplace_back(mesh.vertices.size());
            mesh.vertices.emplace_back(p);
          }
        }
        if (!face_edge_vertices.empty()) {
          if (face_edge_vertices.size() != 2) {
            HERMES_LOG_ERROR("Face intersection not a line {}",
                             face_edge_vertices.size());
          }
          face_vertices.emplace_back(face_edge_vertices[0]);
          face_vertices.emplace_back(face_edge_vertices[1]);
        }
      }
      if (!face_vertices.empty()) {
        // sort edges
        size_t edge_count = face_vertices.size() / 2;
        std::vector<int> visited(edge_count, 0);
        std::vector<size_t> new_face_vertices;
        new_face_vertices.emplace_back(face_vertices[0]);
        new_face_vertices.emplace_back(face_vertices[1]);
        visited[0] = 1;
        size_t next_vertex = face_vertices[1];
        for (int j = 0; j < edge_count; ++j) {
          // find next edge with next_vertex
          for (size_t i = 0; i < edge_count; ++i) {
            if (visited[i])
              continue;
            size_t a = face_vertices[i * 2 + 0];
            size_t b = face_vertices[i * 2 + 1];
            if (a != next_vertex && b != next_vertex)
              continue;
            visited[i] = 1;
            next_vertex = a == next_vertex ? b : a;
            new_face_vertices.emplace_back(next_vertex);
          }
        }
        mesh.faces.emplace_back(new_face_vertices);
      }
    }
    return mesh;
  }

  // *******************************************************************************************************************
  //                                                                                                    PUBLIC FIELDS
  // *******************************************************************************************************************
  std::vector<Patch> patches;
  std::vector<hermes::point3> vertices;
  std::vector<hermes::point3> cell_centers;
  std::vector<hermes::bbox3> cell_bounds;
  std::vector<size_t> owner;
  std::vector<i64> neighbour;
  std::vector<std::vector<size_t>> faces;
  std::vector<std::vector<size_t>> cells;
};

#endif // PSA_ANIM_TOOLS_FOAM2VDB_POLY_MESH_H
