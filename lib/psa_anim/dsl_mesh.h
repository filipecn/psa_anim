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

#ifndef PSA_ANIM_TOOLS_OPENFOAM_PARSER_FA_MESH_H
#define PSA_ANIM_TOOLS_OPENFOAM_PARSER_FA_MESH_H

#include "hermes/logging/logging.h"
#include <boundary_mesh.h>
#include <hermes/geometry/queries.h>
#include <mesh_utils.h>
#include <noise.h>
#include <openfoam_parser.h>
#include <openfoam_poly_mesh.h>
#include <queue>
#include <sim_path.h>
#include <utility>
#include <vector>

/**********************************************************************************************************************/
/*                                                                                                                    */
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
class DSLMesh : public BoundaryMesh {
public:
  // *******************************************************************************************************************
  //                                                                                               DSLMesh::FrontCell
  // *******************************************************************************************************************
  struct FrontCell {
    int cell_id{};
    hermes::point3 position;
    hermes::vec3 velocity;
    real_t height{};
  };
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  ///
  DSLMesh() = default;
  ///
  /// \param poly_mesh openfoam mesh object
  /// \param patch_name openfoam mesh patch name identifier
  explicit DSLMesh(const PolyMesh *poly_mesh,
                   const std::string &patch_name = "terrain")
      : BoundaryMesh(poly_mesh, patch_name) {
    //    set(poly_mesh_, patch_name);
  }
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  ///
  /// \return
  [[nodiscard]] const std::vector<FrontCell> &frontCells() const {
    return front_cells_;
  }
  //                                                                                                       processing
  /// \brief Loads field values from a given frame directory and computes psl
  /// input \note Two extra fields are created in this function: 'front_sdf' and
  /// 'slope' \note This function also calculates the set of front cells \param
  /// path
  void load(const SimPath &sim_path) {
    BoundaryMesh::load(sim_path, {"fa_h", "fa_Us", "fa_W", "fa_divUs"});
    // compute sdf
    scalar_fields_["front_sdf"] = computeFrontDistanceField(front_cells_);
    scalar_fields_["slope"] = computeSlopeField();
  }
  void load(const std::vector<double> &fa_h,
            const std::vector<hermes::vec3> &fa_Us,
            const std::vector<double> &fa_W,
            const std::vector<double> &fa_divUs) {
    scalar_fields_["fa_h"] = fa_h;
    vector_fields_["fa_Us"] = fa_Us;
    scalar_fields_["fa_W"] = fa_W;
    scalar_fields_["fa_divUs"] = fa_divUs;
    scalar_fields_["front_sdf"] = computeFrontDistanceField(front_cells_);
    scalar_fields_["slope"] = computeSlopeField();
  }
  ///
  /// \return
  psa_anim::FaceMesh surfaceMesh(std::vector<double> &front_distance,
                                 std::vector<double> &Ws,
                                 std::vector<double> &mag_Us,
                                 std::vector<double> &u_injs) {

    auto computeUInj = [](double front_distance, const hermes::vec3 &u,
                          double w, const hermes::vec3 &nf,
                          const hermes::point3 &face_center) -> float {
      float psl_front_length_factor = 80;
      float psl_u_entrainment_factor = 1.6;

      // compute front distance decay factor
      auto front_weight = 0.0;

      if (front_distance >= 0 && front_distance < 10000) {
        auto weight = std::exp(front_distance / psl_front_length_factor);
        front_weight = 1.0 / weight;
      }

      float mag_u = u.length();
      float velocity_trigger =
          std::min(1.0, std::max(0.0, 1.0 / (1.0 + std::exp(-mag_u + 10))));

      // w_h is the entrainment velocity
      double entrainment_velocity = 2.0 * std::fabs(w);
      auto u_noise = psa_anim::fractalAllicatorNoise(face_center, 0.5);
      auto a_noise = psa_anim::fractalAllicatorNoise(face_center, 1.0);

      double we = u_noise * mag_u;
      we += entrainment_velocity;
      we *= psl_u_entrainment_factor; // gamma u
      we *= front_weight;
      we *= velocity_trigger;

      // Faces that inject mass must have their boundary conditions updated
      // accordingly
      we = std::max(we, 0.);
      if (we > 0.01) {
        auto normal = (-nf).normalized();
        return ((float)we * normal +
                u * (float)(0.15 * front_weight * velocity_trigger))
            .length();
      }
      return 0.f;
    };

    std::vector<double> dslH_surface_field = scalarField("fa_h");
    std::vector<double> dsl_front_distance_field = scalarField("front_sdf");
    std::vector<double> W_surface_field = scalarField("fa_W");
    std::vector<hermes::vec3> U_surface_field = vectorField("fa_Us");

    // retrieve dsl global index faces (faces with minimal height)
    std::set<size_t> dsl_faces;
    iterateFaces([&](size_t global_face_id, size_t local_face_id) {
      if (dslH_surface_field[local_face_id] < minimal_snow_height)
        return;
      dsl_faces.insert(local_face_id);
    });

    psa_anim::FaceMesh dsl_surface_mesh;
    // get the dsl vertices and remap their indices into the output mesh index
    // global poly mesh index -> output mesh index (vertices vector index)
    std::unordered_map<size_t, size_t> vertex_index_map;
    // here we also register the vertex faces
    // new vertex index -> local face id indices
    std::unordered_map<size_t, std::vector<size_t>> vertex_faces_map;
    // for each dsl surface face register its vertices
    for (size_t face_index : dsl_faces) {
      // get face vertices
      auto face_vertices = faceVertexIndices(face_index);
      // register vertices
      for (auto vertex_index : face_vertices) {
        if (vertex_index_map.find(vertex_index) == vertex_index_map.end()) {
          vertex_index_map[vertex_index] = dsl_surface_mesh.vertices.size();
          dsl_surface_mesh.vertices.emplace_back(vertexAt(vertex_index));
        }
        vertex_faces_map[vertex_index_map[vertex_index]].emplace_back(
            face_index);
      }
    }

    // extrude dsl faces based on neighboring heights
    // also detect boundary vertices
    std::set<size_t> boundary_vertices; // stores new vertex indices
    size_t new_vertex_count = dsl_surface_mesh.vertices.size();
    for (size_t new_vertex_index = 0; new_vertex_index < new_vertex_count;
         ++new_vertex_index) {

      // each vertex is duplicated and moved in the normal's direction
      hermes::vec3 normal;
      real_t height = 0;
      real_t W = 0;
      real_t U = 0;
      real_t front_sdf = 0;
      real_t u_inj_acc = 0;
      size_t front_sdf_count = 0;

      // iterate over neighboring faces
      bool is_boundary = false;

      const auto &neighbour_faces = vertex_faces_map[new_vertex_index];
      for (size_t face_index : neighbour_faces) {
        // accumulate normal and height

        auto face_normal = faceNormal(face_index);
        normal.x += face_normal[0];
        normal.y += face_normal[1];
        normal.z += face_normal[2];
        height += dslH_surface_field[face_index];
        W += W_surface_field[face_index];
        U += U_surface_field[face_index].length();
        if (dsl_front_distance_field[face_index] >= 0 &&
            dsl_front_distance_field[face_index] < 6000) {
          front_sdf += dsl_front_distance_field[face_index];
          front_sdf_count++;

          u_inj_acc += computeUInj(
              dsl_front_distance_field[face_index], U_surface_field[face_index],
              W_surface_field[face_index], face_normal, faceCenter(face_index));
        }

        // if the neighboring face is not a DSL face, then this vertex is a
        // boundary face
        if (!is_boundary) {
          // iterate over each neighbour face and find out
          auto face_neighbours = faceNeighbours(face_index);
          for (auto face_neighbour : face_neighbours)
            if (dsl_faces.find(face_neighbour.face_id) == dsl_faces.end()) {
              is_boundary = true;
              break;
            }
        }
      }
      if (is_boundary)
        boundary_vertices.insert(new_vertex_index);

      HERMES_ASSERT(!neighbour_faces.empty());
      height /= static_cast<real_t>(neighbour_faces.size());
      W /= static_cast<real_t>(neighbour_faces.size());
      U /= static_cast<real_t>(neighbour_faces.size());
      normal.x /= static_cast<real_t>(neighbour_faces.size());
      normal.y /= static_cast<real_t>(neighbour_faces.size());
      normal.z /= static_cast<real_t>(neighbour_faces.size());

      normal.normalize();

      if (front_sdf_count > 0) {
        front_sdf /= static_cast<real_t>(front_sdf_count);
        u_inj_acc /= static_cast<real_t>(front_sdf_count);
      }

      // duplicate vertex (its index MUST be new_vertex_index + new_vertex_count
      HERMES_ASSERT(dsl_surface_mesh.vertices.size() ==
                    new_vertex_index + new_vertex_count);
      // extrude
      dsl_surface_mesh.vertices.emplace_back(
          dsl_surface_mesh.vertices[new_vertex_index] - height * normal);

      front_distance.emplace_back(front_sdf);
      u_injs.emplace_back(u_inj_acc);
      mag_Us.emplace_back(U);
      Ws.emplace_back(W);
    }

    // generate new mesh faces
    for (size_t face_index : dsl_faces) {
      auto face_vertices = faceVertexIndices(face_index);
      // convert vertex indices to new indices
      for (auto &v : face_vertices) {
        HERMES_ASSERT(vertex_index_map.count(v))
        v = vertex_index_map[v];
      }
      // terrain face
      dsl_surface_mesh.faces.emplace_back(face_vertices);
      // surface face
      {
        std::vector<size_t> vertex_indices;
        for (auto vertex_index : face_vertices)
          vertex_indices.emplace_back(vertex_index + new_vertex_count);
        dsl_surface_mesh.faces.emplace_back(vertex_indices);
      }
      // boundary faces
      {
        auto isBoundary = [&](size_t face_vertex) -> bool {
          return boundary_vertices.find(face_vertex) != boundary_vertices.end();
        };
        // detect boundary edges
        for (size_t i = 0; i < face_vertices.size(); ++i) {
          size_t va = face_vertices[i];
          size_t vb = face_vertices[(i + 1) & face_vertices.size()];
          if (isBoundary(va) && isBoundary(vb)) {
            // create boundary face
            std::vector<size_t> vertex_indices = {vb, va, va + new_vertex_count,
                                                  vb + new_vertex_count};
            //            dsl_surface_mesh.faces.emplace_back(vertex_indices);
          }
        }
      }
    }

    return std::move(dsl_surface_mesh);
  }
  // *******************************************************************************************************************
  //                                                                                                    PUBLIC FIELDS
  // *******************************************************************************************************************
  double minimal_snow_height = 0.02;
  double minimal_velocity = 0.01;
  double minimal_snow_height_ratio = 1.5;

private:
  // *******************************************************************************************************************
  //                                                                                                        FRONT SDF
  // *******************************************************************************************************************

  std::vector<double>
  computeFrontDistanceField2(std::vector<FrontCell> &front_face_list) {
    // init face sdf distance
    std::vector<double> sdf_field_ref(poly_mesh_->patches[patch_id_].size, -1);
    // read dsl data
    HERMES_ASSERT(poly_mesh_->patch(patch_name_).size);
    HERMES_ASSERT(scalar_fields_.count("fa_h"));
    HERMES_ASSERT(vector_fields_.count("fa_Us"));
    const std::vector<hermes::vec3> &dslU_surface_field =
        vector_fields_["fa_Us"];
    std::vector<double> dslH_surface_field = scalar_fields_["fa_h"];
    HERMES_ASSERT(poly_mesh_->patch(patch_name_).size ==
                  dslH_surface_field.size());
    HERMES_ASSERT(poly_mesh_->patch(patch_name_).size ==
                  dslU_surface_field.size());

    // detect faces that contain the iso-line
    // lead edge faces (faces containing the front iso-line) are faces on the
    // very front edge of the flow, adjacent to empty cells in their velocity
    // direction
    // this set stores global face ids
    std::set<size_t> candidate_front_faces;
    for (size_t local_face_id = 0;
         local_face_id < poly_mesh_->patches[patch_id_].size; ++local_face_id) {
      auto global_face_id =
          local_face_id + poly_mesh_->patches[patch_id_].start;
      // get face info
      const auto face_h = dslH_surface_field[local_face_id];
      const auto face_U = dslU_surface_field[local_face_id];
      const auto face_center = poly_mesh_->faceCenter(global_face_id);

      // filter dsl faces
      if (face_h < minimal_snow_height || face_U.length() < minimal_velocity)
        continue;
      // retrieve neighbouring faces
      auto neighbours = globalFaceNeighbours(global_face_id);
      HERMES_ASSERT(!neighbours.empty());
      // iterate over neighbor faces
      for (auto neighbour : neighbours) {
        // Method 1
        // Consider large gradients along velocity direction DOWNWARDS
        const auto neighbour_center = poly_mesh_->faceCenter(neighbour.face_id);
        // check if neighbour is towards velocity direction

        // here we see if this neighbour cell is ahead of the center cell
        auto edge_a = vertexAt(neighbour.edge_a).xy();
        auto edge_b = vertexAt(neighbour.edge_b).xy();
        hermes::Line2 edge_line(edge_a, edge_b - edge_a);
        if (hermes::GeometricQueries::intersect(
                edge_line, hermes::ray2(face_center.xy(), face_U.xy())) &&
            neighbour_center.z <= face_center.z) {
          const auto neighbour_local_face_id =
              neighbour.face_id - poly_mesh_->patches[patch_id_].start;
          const auto neighbour_h = dslH_surface_field[neighbour_local_face_id];
          // check h gradient
          if (neighbour_h < minimal_snow_height ||
              (face_h / neighbour_h > minimal_snow_height_ratio &&
               face_U.length() >
                   dslU_surface_field[neighbour_local_face_id].length())) {
            candidate_front_faces.insert(global_face_id);
            break;
          }
        }
      }
    }

    // candidate_front_faces potentially has too much faces
    // here we filter the candidate faces to a list of actual front faces
    std::queue<int> lead_edge_faces;
    front_face_list.clear();
    // real front faces will be the faces with no candidates ahead
    for (const auto &candidate_front_face : candidate_front_faces) {
      // indices
      auto global_face_id = candidate_front_face;
      HERMES_ASSERT(global_face_id >= poly_mesh_->patches[patch_id_].start);
      auto local_face_id =
          global_face_id - poly_mesh_->patches[patch_id_].start;
      // data
      const auto face_U = dslU_surface_field[local_face_id];
      const auto face_center = poly_mesh_->faceCenter(global_face_id);

      // get neighbors along the velocity direction
      auto neighbours = globalFaceNeighbours(global_face_id);
      HERMES_ASSERT(!neighbours.empty());

      // iterate over neighbor faces
      bool is_leading_edge = true;
      for (auto neighbour : neighbours) {
        // Method 1
        // Consider large gradients along velocity direction DOWNWARDS
        const auto neighbour_center = poly_mesh_->faceCenter(neighbour.face_id);
        // check if neighbour is towards velocity direction

        auto edge_a = vertexAt(neighbour.edge_a).xy();
        auto edge_b = vertexAt(neighbour.edge_b).xy();
        hermes::Line2 edge_line(edge_a, edge_b - edge_a);
        if (!hermes::GeometricQueries::intersect(
                edge_line, hermes::ray2(face_center.xy(), face_U.xy())) &&
            candidate_front_faces.find(neighbour.face_id) !=
                candidate_front_faces.end()) {
          is_leading_edge = false;
          break;
        }
      }

      is_leading_edge = true;
      if (is_leading_edge) {
        sdf_field_ref[local_face_id] = 0;
        lead_edge_faces.push(local_face_id);
        front_face_list.push_back(
            {.cell_id = (int)local_face_id,
             .position = poly_mesh_->faceCenter(
                 local_face_id + poly_mesh_->patches[patch_id_].start),
             .velocity = dslU_surface_field[local_face_id],
             .height = static_cast<real_t>(dslH_surface_field[local_face_id])});
      }
    }

    // compute sdf from isoline
    while (!lead_edge_faces.empty()) {
      auto local_face_id = lead_edge_faces.front();
      lead_edge_faces.pop();
      auto global_face_id =
          local_face_id + poly_mesh_->patches[patch_id_].start;
      const auto face_U = dslU_surface_field[local_face_id];
      const auto face_center = poly_mesh_->faceCenter(global_face_id);
      // consider only neighbors on the opposite direction of the velocity
      // direction retrieve neighbouring faces iterate over boundary neighbor
      // faces
      auto neighbours = globalFaceNeighbours(global_face_id);
      for (auto neighbour : neighbours) {
        const auto neighbour_center = poly_mesh_->faceCenter(neighbour.face_id);
        auto edge_a = vertexAt(neighbour.edge_a).xy();
        auto edge_b = vertexAt(neighbour.edge_b).xy();
        hermes::Line2 edge_line(edge_a, edge_b - edge_a);
        if (hermes::GeometricQueries::intersect(
                edge_line, hermes::ray2(face_center.xy(), -face_U.xy()))) {
          const auto neighbour_local_face_id =
              neighbour.face_id - poly_mesh_->patches[patch_id_].start;
          const auto neighbour_h = dslH_surface_field[neighbour_local_face_id];
          // discard no dsl faces
          if (neighbour_h < minimal_snow_height)
            continue;
          // push face only if it can be updated
          auto dist = (neighbour_center - face_center).length();
          // propagate only to greater values
          if (sdf_field_ref[neighbour_local_face_id] >= 0 &&
              sdf_field_ref[neighbour_local_face_id] <=
                  sdf_field_ref[local_face_id] + dist)
            continue;
          // little hack due to -1 values for outside dsl faces
          auto fixed_init_value = sdf_field_ref[neighbour_local_face_id] < 0
                                      ? std::numeric_limits<double>::max()
                                      : sdf_field_ref[neighbour_local_face_id];
          if (fixed_init_value > sdf_field_ref[local_face_id] + dist) {
            sdf_field_ref[neighbour_local_face_id] =
                sdf_field_ref[local_face_id] + dist;
            lead_edge_faces.push(neighbour_local_face_id);
          }
        }
      }
    }
    // convert -1 values into max values
    for (auto &value : sdf_field_ref)
      if (value < 0)
        value = hermes::Numbers::greatest_f32();
    return sdf_field_ref;
  }

  /// Computes front distance field by front propagation
  /// \param front_face_list output list of leading edge faces
  /// \return A vector containing the distance from the front for every cell in
  /// the patch
  std::vector<double>
  computeFrontDistanceField(std::vector<FrontCell> &front_face_list) {
    // init face sdf distance
    std::vector<double> sdf_field_ref(poly_mesh_->patches[patch_id_].size, -1);
    // read dsl data
    HERMES_ASSERT(poly_mesh_->patch(patch_name_).size);
    HERMES_ASSERT(scalar_fields_.count("fa_h"));
    HERMES_ASSERT(vector_fields_.count("fa_Us"));
    const std::vector<hermes::vec3> &dslU_surface_field =
        vector_fields_["fa_Us"];
    std::vector<double> dslH_surface_field = scalar_fields_["fa_h"];
    HERMES_ASSERT(poly_mesh_->patch(patch_name_).size ==
                  dslH_surface_field.size());
    HERMES_ASSERT(poly_mesh_->patch(patch_name_).size ==
                  dslU_surface_field.size());

    // detect faces that contain the iso-line
    // lead edge faces (faces containing the front iso-line) are faces on the
    // very front edge of the flow, adjacent to empty cells in their velocity
    // direction
    // this set stores global face ids
    std::set<size_t> candidate_front_faces;
    for (size_t local_face_id = 0;
         local_face_id < poly_mesh_->patches[patch_id_].size; ++local_face_id) {
      auto global_face_id =
          local_face_id + poly_mesh_->patches[patch_id_].start;
      // get face info
      const auto face_h = dslH_surface_field[local_face_id];
      const auto face_U = dslU_surface_field[local_face_id];
      const auto face_center = poly_mesh_->faceCenter(global_face_id);

      // filter dsl faces
      if (face_h < minimal_snow_height || face_U.length() < minimal_velocity)
        continue;
      // retrieve neighbouring faces
      auto neighbours = globalFaceNeighbours(global_face_id);
      HERMES_ASSERT(!neighbours.empty());
      // iterate over neighbor faces
      for (auto neighbour : neighbours) {
        // Method 1
        // Consider large gradients along velocity direction DOWNWARDS
        const auto neighbour_center = poly_mesh_->faceCenter(neighbour.face_id);
        // check if neighbour is towards velocity direction

        // here we see if this neighbour cell is ahead of the center cell
        auto edge_a = vertexAt(neighbour.edge_a).xy();
        auto edge_b = vertexAt(neighbour.edge_b).xy();
        hermes::Line2 edge_line(edge_a, edge_b - edge_a);
        if (hermes::GeometricQueries::intersect(
                edge_line, hermes::ray2(face_center.xy(), face_U.xy()))) {
          const auto neighbour_local_face_id =
              neighbour.face_id - poly_mesh_->patches[patch_id_].start;
          const auto neighbour_h = dslH_surface_field[neighbour_local_face_id];
          // check h gradient
          if (neighbour_h < minimal_snow_height ||
              (face_h / neighbour_h > minimal_snow_height_ratio &&
               face_U.length() >
                   dslU_surface_field[neighbour_local_face_id].length())) {
            candidate_front_faces.insert(global_face_id);
            break;
          }
        }
      }
    }

    // candidate_front_faces potentially has too much faces
    // here we filter the candidate faces to a list of actual front faces
    std::queue<int> lead_edge_faces;
    front_face_list.clear();
    // real front faces will be the faces with no candidates ahead
    for (const auto &candidate_front_face : candidate_front_faces) {
      // indices
      auto global_face_id = candidate_front_face;
      HERMES_ASSERT(global_face_id >= poly_mesh_->patches[patch_id_].start);
      auto local_face_id =
          global_face_id - poly_mesh_->patches[patch_id_].start;
      // data
      const auto face_U = dslU_surface_field[local_face_id];
      const auto face_center = poly_mesh_->faceCenter(global_face_id);

      // get neighbors along the velocity direction
      auto neighbours = globalFaceNeighbours(global_face_id);
      HERMES_ASSERT(!neighbours.empty());

      // iterate over neighbor faces
      bool is_leading_edge = true;
      for (auto neighbour : neighbours) {
        // Method 1
        // Consider large gradients along velocity direction DOWNWARDS
        const auto neighbour_center = poly_mesh_->faceCenter(neighbour.face_id);
        // check if neighbour is towards velocity direction

        auto edge_a = vertexAt(neighbour.edge_a).xy();
        auto edge_b = vertexAt(neighbour.edge_b).xy();
        hermes::Line2 edge_line(edge_a, edge_b - edge_a);
        if (!hermes::GeometricQueries::intersect(
                edge_line, hermes::ray2(face_center.xy(), face_U.xy())) &&
            candidate_front_faces.find(neighbour.face_id) !=
                candidate_front_faces.end()) {
          is_leading_edge = false;
          break;
        }
      }

      is_leading_edge = true;
      if (is_leading_edge) {
        sdf_field_ref[local_face_id] = 0;
        lead_edge_faces.push(local_face_id);
        front_face_list.push_back(
            {.cell_id = (int)local_face_id,
             .position = poly_mesh_->faceCenter(
                 local_face_id + poly_mesh_->patches[patch_id_].start),
             .velocity = dslU_surface_field[local_face_id],
             .height = static_cast<real_t>(dslH_surface_field[local_face_id])});
      }
    }

    // compute sdf from isoline
    while (!lead_edge_faces.empty()) {
      auto local_face_id = lead_edge_faces.front();
      lead_edge_faces.pop();
      auto global_face_id =
          local_face_id + poly_mesh_->patches[patch_id_].start;
      const auto face_U = dslU_surface_field[local_face_id];
      const auto face_center = poly_mesh_->faceCenter(global_face_id);
      // consider only neighbors on the opposite direction of the velocity
      // direction retrieve neighbouring faces iterate over boundary neighbor
      // faces
      auto neighbours = globalFaceNeighbours(global_face_id);
      for (auto neighbour : neighbours) {
        const auto neighbour_center = poly_mesh_->faceCenter(neighbour.face_id);
        auto edge_a = vertexAt(neighbour.edge_a).xy();
        auto edge_b = vertexAt(neighbour.edge_b).xy();
        hermes::Line2 edge_line(edge_a, edge_b - edge_a);
        if (hermes::GeometricQueries::intersect(
                edge_line, hermes::ray2(face_center.xy(), -face_U.xy()))) {
          const auto neighbour_local_face_id =
              neighbour.face_id - poly_mesh_->patches[patch_id_].start;
          const auto neighbour_h = dslH_surface_field[neighbour_local_face_id];
          // discard no dsl faces
          if (neighbour_h < minimal_snow_height)
            continue;
          // push face only if it can be updated
          auto dist = (neighbour_center - face_center).length();
          // propagate only to greater values
          if (sdf_field_ref[neighbour_local_face_id] >= 0 &&
              sdf_field_ref[neighbour_local_face_id] <=
                  sdf_field_ref[local_face_id] + dist)
            continue;
          // little hack due to -1 values for outside dsl faces
          auto fixed_init_value = sdf_field_ref[neighbour_local_face_id] < 0
                                      ? std::numeric_limits<double>::max()
                                      : sdf_field_ref[neighbour_local_face_id];
          if (fixed_init_value > sdf_field_ref[local_face_id] + dist) {
            sdf_field_ref[neighbour_local_face_id] =
                sdf_field_ref[local_face_id] + dist;
            lead_edge_faces.push(neighbour_local_face_id);
          }
        }
      }
    }
    // convert -1 values into max values
    for (auto &value : sdf_field_ref)
      if (value < 0)
        value = hermes::Numbers::greatest_f32();
    return sdf_field_ref;
  }
  // *******************************************************************************************************************
  //                                                                                                            SLOPE
  // *******************************************************************************************************************
  /// Computes the terrain slope for each cell in the selected poly mesh patch
  std::vector<double> computeSlopeField() {
    HERMES_ASSERT(poly_mesh_->patch(patch_name_).size);
    std::vector<double> values;
    size_t n = poly_mesh_->patch(patch_name_).size;
    for (size_t i = 0; i < n; ++i) {
      // compute slope
      auto face_normal = hermes::normalize(
          -poly_mesh_->faceNormal(poly_mesh_->patches[patch_id_].start + i));
      values.emplace_back(std::acos(face_normal.z));
    }
    return values;
  }
  // *******************************************************************************************************************
  //                                                                                                           FIELDS
  // *******************************************************************************************************************
  //                                                                                                             data
  std::vector<FrontCell> front_cells_;
  std::vector<hermes::point3> cell_centers_;
};

#endif // PSA_ANIM_TOOLS_OPENFOAM_PARSER_FA_MESH_H
