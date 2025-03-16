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
///\file fa_mesh.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-11-03
///
///\brief

#ifndef PSA_ANIM_BOUNDARY_MESH_H
#define PSA_ANIM_BOUNDARY_MESH_H

#include <hermes/geometry/queries.h>
#include <openfoam_parser.h>
#include <openfoam_poly_mesh.h>
#include <sim_path.h>
#include <utility>
#include <vector>

// *********************************************************************************************************************
//                                                                                                       BoundaryMesh
// *********************************************************************************************************************
class BoundaryMesh {
public:
  // *******************************************************************************************************************
  //                                                                                          BoundaryMesh::Neighbour
  // *******************************************************************************************************************
  struct Neighbour {
    size_t face_id;
    size_t edge_a;
    size_t edge_b;
  };
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  ///
  BoundaryMesh() = default;
  ///
  /// \param poly_mesh openfoam mesh object
  /// \param patch_name openfoam mesh patch name identifier
  BoundaryMesh(const PolyMesh *poly_mesh, const std::string &patch_name)
      : poly_mesh_(poly_mesh) {
    set(poly_mesh_, patch_name);
  }
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  ///
  /// \return
  [[nodiscard]] bool good() const {
    return poly_mesh_ && patch_id_ >= 0 &&
           patch_id_ < poly_mesh_->patches.size();
  }
  //                                                                                                           access
  void iterateFaces(std::function<void(size_t, size_t)> callback) {
    HERMES_LOG_AND_RETURN_IF_NOT(good(), "Bad patch.");
    auto start = poly_mesh_->patches[patch_id_].start;
    auto end = start + poly_mesh_->patches[patch_id_].size;
    for (size_t global_face_id = start; global_face_id < end; ++global_face_id)
      callback(global_face_id, global_face_id - start);
  }
  ///
  /// \return
  [[nodiscard]] std::string patchName() const {
    HERMES_ASSERT(good());
    return poly_mesh_->patches[patch_id_].name;
  }
  ///
  /// \param global_face_id
  /// \return
  [[nodiscard]] size_t localFaceId(size_t global_face_id) const {
    HERMES_ASSERT(good());
    HERMES_ASSERT(global_face_id >= poly_mesh_->patches[patch_id_].start);
    HERMES_ASSERT(global_face_id < poly_mesh_->patches[patch_id_].start +
                                       poly_mesh_->patches[patch_id_].size);
    return global_face_id - poly_mesh_->patches[patch_id_].start;
  }
  ///
  /// \param local_face_id
  /// \return
  [[nodiscard]] size_t globalFaceId(size_t local_face_id) const {
    HERMES_ASSERT(good());
    HERMES_ASSERT(local_face_id < size());
    return local_face_id + poly_mesh_->patches[patch_id_].start;
  }
  ///
  /// \return
  [[nodiscard]] size_t size() const {
    if (!good())
      return 0;
    return poly_mesh_->patches[patch_id_].size;
  }
  ///
  /// \param local_face_id
  /// \return
  [[nodiscard]] size_t faceSize(size_t local_face_id) const {
    size_t global_face_id = globalFaceId(local_face_id);
    return poly_mesh_->faceSize(global_face_id);
  }
  ///
  /// \param local_face_id
  /// \return
  [[nodiscard]] hermes::vec3 faceNormal(size_t local_face_id) const {
    return poly_mesh_->faceNormal(globalFaceId(local_face_id));
  }

  ///
  /// \param local_face_id
  /// \return local face ids
  [[nodiscard]] std::vector<Neighbour>
  faceNeighbours(size_t local_face_id) const {
    HERMES_ASSERT(local_face_id < poly_mesh_->patches[patch_id_].size);
    auto global_face_id = local_face_id + poly_mesh_->patches[patch_id_].start;
    std::vector<Neighbour> neighbours;
    auto ns = globalFaceNeighbours(global_face_id);
    for (const auto &n : ns)
      neighbours.push_back({.face_id = localFaceId(n.face_id),
                            .edge_a = n.edge_a,
                            .edge_b = n.edge_b});
    return neighbours;
  }
  ///
  /// \return
  [[nodiscard]] const std::vector<hermes::point3> &faceCenters() const {
    return cell_centers_;
  }

  [[nodiscard]] hermes::point3 faceCenter(size_t local_face_id) const {
    return poly_mesh_->faceCenter(globalFaceId(local_face_id));
  }
  ///
  /// \param local_face_id
  /// \param face_vertex_index
  /// \return
  [[nodiscard]] const hermes::point3 &
  faceVertex(size_t local_face_id, size_t face_vertex_index) const {
    size_t global_face_id = globalFaceId(local_face_id);
    return poly_mesh_
        ->vertices[poly_mesh_->faces[global_face_id][face_vertex_index]];
  }

  [[nodiscard]] const std::vector<size_t> &
  faceVertexIndices(size_t local_face_id) const {
    return poly_mesh_->faces[globalFaceId(local_face_id)];
  }
  ///
  /// \param global_vertex_index
  /// \return
  [[nodiscard]] const hermes::point3 &
  vertexAt(size_t global_vertex_index) const {
    HERMES_ASSERT(global_vertex_index < poly_mesh_->vertices.size());
    return poly_mesh_->vertices[global_vertex_index];
  }
  //                                                                                                           fields
  ///
  /// \param field_name
  /// \return
  [[nodiscard]] bool hasField(const std::string &field_name) const {
    return scalar_fields_.count(field_name) || vector_fields_.count(field_name);
  }
  ///
  /// \param field_name
  /// \return
  [[nodiscard]] const std::vector<double> &
  scalarField(const std::string &field_name) const {
    static std::vector<double> dummy;
    auto it = scalar_fields_.find(field_name);
    if (it != scalar_fields_.end())
      return it->second;
    return dummy;
  }
  ///
  /// \param field_name
  /// \param local_face_id
  /// \return
  [[nodiscard]] double scalarAt(const std::string &field_name,
                                size_t local_face_id) const {
    auto it = scalar_fields_.find(field_name);
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(it != scalar_fields_.end(), 0,
                                       "field not found");
    return it->second[local_face_id];
  }
  ///
  /// \param field_name
  /// \param local_face_id
  /// \return
  [[nodiscard]] hermes::vec3 vectorAt(const std::string &field_name,
                                      size_t local_face_id) const {
    auto it = vector_fields_.find(field_name);
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(it != vector_fields_.end(), {},
                                       "field not found");
    return it->second[local_face_id];
  }
  ///
  /// \param field_name
  /// \return
  [[nodiscard]] double scalarFieldMin(const std::string &field_name) const {
    auto it = scalar_fields_.find(field_name);
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(it != scalar_fields_.end(), 0,
                                       "invalid scalar field name");
    const auto &values = it->second;
    f64 min_value = hermes::Numbers::greatest<f64>();
    for (auto value : values) {
      min_value = std::min(min_value, value);
    }
    return min_value;
  }
  ///
  /// \param field_name
  /// \return
  [[nodiscard]] double scalarFieldMax(const std::string &field_name) const {
    auto it = scalar_fields_.find(field_name);
    HERMES_LOG_AND_RETURN_VALUE_IF_NOT(it != scalar_fields_.end(), 0,
                                       "invalid scalar field name");
    const auto &values = it->second;
    f64 max_value = hermes::Numbers::lowest_f64();
    for (auto value : values) {
      max_value = std::max(max_value, value);
    }
    return max_value;
  }
  ///
  /// \param field_name
  /// \return
  [[nodiscard]] const std::vector<hermes::vec3> &
  vectorField(const std::string &field_name) const {
    static std::vector<hermes::vec3> dummy;
    auto it = vector_fields_.find(field_name);
    if (it != vector_fields_.end())
      return it->second;
    return dummy;
  }
  ///
  /// \return
  [[nodiscard]] std::vector<std::string> scalarFieldNames() const {
    std::vector<std::string> names;
    for (const auto &field : scalar_fields_)
      names.emplace_back(field.first);
    return names;
  }
  ///
  /// \return
  [[nodiscard]] std::vector<std::string> vectorFieldNames() const {
    std::vector<std::string> names;
    for (const auto &field : vector_fields_)
      names.emplace_back(field.first);
    return names;
  }
  ///
  /// \param name
  /// \param values
  void pushScalarField(const std::string &name,
                       const std::vector<double> &values) {
    // HERMES_LOG("pushing field {} with {} values", name, values.size());
    scalar_fields_[name] = values;
  }
  ///
  /// \param name
  /// \param values
  void pushVectorField(const std::string &name,
                       const std::vector<hermes::vec3> &values) {
    // HERMES_LOG("pushing field {} with {} values", name, values.size());
    vector_fields_[name] = values;
  }
  //                                                                                                     construction
  ///
  /// \param poly_mesh
  /// \param patch_name
  void set(const PolyMesh *poly_mesh, const std::string &patch_name) {
    poly_mesh_ = poly_mesh;
    patch_name_ = patch_name;
    HERMES_LOG_AND_RETURN_IF_NOT(poly_mesh->containsPatch(patch_name),
                                 "patch not found");
    patch_id_ = poly_mesh_->patchId(patch_name);
    buildTopology();
  }
  //                                                                                                       processing
  /// \brief Loads field values from a given frame directory and computes psl
  /// input
  /// \note Two extra fields are created in this function: 'front_sdf' and
  /// 'slope'
  /// \note This function also calculates the set of front cells
  /// \param path
  void load(const SimPath &sim_path,
            const std::vector<std::string> &field_names_filter = {}) {
    auto path = sim_path.currentFramePath();
    HERMES_LOG_AND_RETURN_IF_NOT(good(), "patch not set");
    // HERMES_LOG("Loading Boundary data from {}", path);
    HERMES_LOG_AND_RETURN_IF_NOT(path.exists() && path.isDirectory(),
                                 "Invalid frame path");
    // cleanup previous data
    scalar_fields_.clear();
    vector_fields_.clear();
    std::vector<std::string> scalar_field_names, vector_field_names;
    sim_path.listFields(scalar_field_names, vector_field_names);
    // load frame's data
    if (field_names_filter.empty()) {
      // load all fields
      for (const auto &name : scalar_field_names)
        readField(path, name);
      for (const auto &name : vector_field_names)
        readField(path, name);
    } else {
      for (const auto &s : field_names_filter)
        readField(path, s);
    }
  }
  // *******************************************************************************************************************
  //                                                                                                    PUBLIC FIELDS
  // *******************************************************************************************************************
protected:
  /// Reconstructs OpenFOAM's poly mesh topology by populating a map of edges
  /// Edges represent the connection of two neighboring faces and form a graph
  /// from the mesh faces (the dual mesh) Thus, the map of edges maps a pair of
  /// primal vertices (corresponding to a primal cell face) to a pair of
  /// vertices in the dual structure, which corresponds to a pair of cell ids
  void buildTopology() {
    HERMES_LOG("Building poly mesh patch topology...");
    edges_.clear();
    const size_t face_start = poly_mesh_->patches[patch_id_].start;
    const size_t face_end = face_start + poly_mesh_->patches[patch_id_].size;
    // Iterate over all cells and produce an dual edge for each cell face
    // connecting both cells
    cell_centers_.clear();
    for (size_t face_id = face_start; face_id < face_end; ++face_id) {
      cell_centers_.emplace_back(poly_mesh_->faceCenter(face_id));
      auto edge_keys = faceEdgeKeys(face_id);
      for (const auto k : edge_keys) {
        if (edges_.count(k))
          // this is the second face for this edge
          edges_[k].second = face_id;
        else
          edges_[k].first = edges_[k].second = face_id;
      }
    }
  }
  /// Reads field data for this patch
  /// \note the field type must be volScalarField or volVectorField
  /// \param input file field input path
  /// \param output_path output file path
  /// \param field_name
  void readField(const hermes::Path &input_directory,
                 const std::string &field_name) {
    auto input_file = input_directory / field_name;
    // HERMES_LOG("loading field from {}", input_file);
    HERMES_LOG_AND_RETURN_IF_NOT(input_file.isFile(),
                                 "could not open input file.")

    OpenFoamDict dict(input_file);
    HERMES_ASSERT(dict.nodes().count("FoamFile"));
    auto type = dict["FoamFile"]["class"].value;
    // HERMES_LOG_VARIABLE(type);
    // HERMES_ASSERT(dict.nodes().count("boundaryField"));
    // HERMES_ASSERT(dict["boundaryField"].fields.count(patch_name_));
    if (dict["boundaryField"][patch_name_].fields.count("value") == 0) {
      HERMES_LOG_WARNING("no value found for field {}", field_name);
      return;
    }
    if (type == "volScalarField") {
      scalar_fields_[field_name] = OpenFoamDict::parseValuesFrom<double>(
          dict["boundaryField"][patch_name_]["value"].value);
      // HERMES_LOG("{} values read", scalar_fields_[field_name].size());
    } else if (type == "volVectorField") {
      vector_fields_[field_name] = OpenFoamDict::parseValuesFrom<hermes::vec3>(
          dict["boundaryField"][patch_name_]["value"].value);
      // HERMES_LOG("{} values read", vector_fields_[field_name].size());
    } else
      HERMES_LOG_ERROR("invalid field type");
  }
  // *******************************************************************************************************************
  //                                                                                                AUXILIARY METHODS
  // *******************************************************************************************************************
  /// Computes the edge map key for a given face
  /// \param global_face_id OpenFOAM's global index for the given face
  /// \return the primal (vertex,vertex) pair corresponding to a primal face
  [[nodiscard]] std::vector<std::pair<size_t, size_t>>
  faceEdgeKeys(size_t global_face_id) const {
    std::vector<std::pair<size_t, size_t>> edge_keys;
    auto key = [](size_t a, size_t b) {
      return std::make_pair(std::min(a, b), std::max(a, b));
    };
    // for each face, register its edges
    const size_t n = poly_mesh_->faces[global_face_id].size();
    for (size_t i = 0; i < n; ++i) {
      size_t a = poly_mesh_->faces[global_face_id][i];
      size_t b = poly_mesh_->faces[global_face_id][(i + 1) % n];
      edge_keys.emplace_back(key(a, b));
    }
    return edge_keys;
  }
  /// Retrieves the list of edges for a given cell
  /// \param global_face_id global poly mesh cell id
  /// \param neighbour_edges
  /// \return list of neighboring global face ids
  [[nodiscard]] std::vector<Neighbour>
  globalFaceNeighbours(size_t global_face_id) const {
    std::vector<Neighbour> neighbours;
    auto neighbourFace = [&](const std::pair<size_t, size_t> &key) {
      return global_face_id == key.first ? key.second : key.first;
    };
    auto keys = faceEdgeKeys(global_face_id);
    for (const auto &key : keys) {
      auto it = edges_.find(key);
      if (it != edges_.end()) {
        auto n = neighbourFace(it->second);
        if (n != global_face_id) {
          HERMES_ASSERT(it->first.first < poly_mesh_->vertices.size());
          HERMES_ASSERT(it->first.second < poly_mesh_->vertices.size());
          neighbours.push_back({.face_id = n,
                                .edge_a = it->first.first,
                                .edge_b = it->first.second});
        }
      } else {
        HERMES_ASSERT(edges_.count(key));
        exit(-1);
      }
    }
    return neighbours;
  }
  // *******************************************************************************************************************
  //                                                                                                           FIELDS
  // *******************************************************************************************************************
  //                                                                                                        poly mesh
  const PolyMesh *poly_mesh_{nullptr};
  std::string patch_name_;
  int patch_id_{-1};
  //                                                                                                             data
  std::map<std::string, std::vector<double>> scalar_fields_;
  std::map<std::string, std::vector<hermes::vec3>> vector_fields_;
  std::vector<hermes::point3> cell_centers_;
  //                                                                                                         topology
  /// Auxiliary structure to represent the dual mesh (graph) of the terrain mesh
  /// Here, a pair of vertices from the primal mesh correspond to a primal face
  /// which in turn corresponds to an edge in the dual structure, which connects
  /// two faces
  // vertex, vertex -> global face, global face
  std::map<std::pair<size_t, size_t>, std::pair<size_t, size_t>> edges_;
};

#endif // PSA_ANIM_TOOLS_OPENFOAM_PARSER_FA_MESH_H
