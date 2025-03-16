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
///\file block_mesh.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-11-30
///
///\brief

#include <block_mesh.h>

namespace psa_anim {

BlockMeshDescription
blockMesh(const std::vector<hermes::point3> &vertices,
          const std::vector<std::vector<size_t>> &faces,
          std::map<std::string, std::vector<hermes::vec3>> boundaries,
          real_t height, BlockMeshConfig config, hermes::vec2 slope_direction,
          hermes::point2 slope_point, bool top_is_inclined) {
  HERMES_LOG("Generating blockMesh dict");
  HERMES_LOG("==================================================");
  HERMES_LOG_VARIABLE(height);
  HERMES_LOG_VARIABLE((int)config.face_matching);
  HERMES_LOG_VARIABLE(config.grading);
  HERMES_LOG_VARIABLE(config.resolution);
  HERMES_LOG_VARIABLE(slope_direction);
  HERMES_LOG_VARIABLE(slope_point);
  HERMES_LOG_VARIABLE((int)top_is_inclined);
  HERMES_LOG("==================================================");

  // return
  BlockMeshDescription block_mesh;
  FaceMesh bad_faces;

  // check boundary input
  HERMES_ASSERT(boundaries.count("bottom") || boundaries.count("terrain"));
  HERMES_ASSERT(boundaries.count("top"));
  std::string top_patch_name = "top";
  std::string bottom_patch_name =
      boundaries.count("bottom") ? "bottom" : "terrain";

  // check if input mesh is quad
  for (const auto &face : faces)
    HERMES_ASSERT(face.size() == 4);

  // extrude, but no need to generate boundary faces
  auto input_mesh = extrudeZ({vertices, faces}, height, slope_direction,
                             slope_point, top_is_inclined, false);

  bad_faces.vertices = input_mesh.vertices;

  // after extrusion, the first half of face should point downwards and the
  // second half upwards
  HERMES_ASSERT(input_mesh.faces.size() == 2 * faces.size());
  for (size_t i = 0; i < faces.size(); ++i) {
    auto down_normal = faceNormal(input_mesh.vertices, input_mesh.faces[i]);
    auto up_normal =
        faceNormal(input_mesh.vertices, input_mesh.faces[i + faces.size()]);
    auto down_alignment =
        directionAlignment(down_normal, hermes::vec3(0, 0, -1));
    auto inverse_down_alignment =
        directionAlignment(down_normal, hermes::vec3(0, 0, 1));
    auto up_alignment = directionAlignment(up_normal, hermes::vec3(0, 0, 1));
    auto inverse_up_alignment =
        directionAlignment(up_normal, hermes::vec3(0, 0, -1));
    HERMES_ASSERT(down_alignment > inverse_down_alignment);
    HERMES_ASSERT(up_alignment > inverse_up_alignment);
  }

  // before going any further, lets compute boundary data
  auto edge_map = computeEdgeTopology(faces);
  // wall boundary names
  std::vector<std::string> wall_patch_names;
  for (const auto &boundary : boundaries)
    if (boundary.first != bottom_patch_name && boundary.first != top_patch_name)
      wall_patch_names.emplace_back(boundary.first);
  HERMES_LOG_ARRAY(wall_patch_names);

  // this function matches a given face with the boundary label
  auto labelDirection =
      [&](const hermes::vec3 &d,
          const std::vector<std::string> &filter) -> std::string {
    std::string largest_id;
    float largest_align = -1;
    for (const auto &name : filter) {
      HERMES_ASSERT(boundaries.count(name));
      for (const auto &direction : boundaries[name]) {
        auto a = directionAlignment(direction, d);
        if (a > largest_align) {
          largest_align = a;
          largest_id = name;
        }
      }
    }
    if (largest_align <= 0)
      return "";
    return largest_id;
  };

  // create block mesh
  for (const auto &vertex : input_mesh.vertices)
    block_mesh.pushVertex(vertex);

  // blocks follow the vertex sequence
  //                 7 _______  6
  //                /|        /
  //              4 _______5/  |
  //              |  |    |   |
  //  x3   x2     |   3 --|-- / 2
  //  | /         | /     | /
  //   -- x1      0  ------ 1

  // notice that the whole ordering is based on the origin of the local
  // coordinate system which is based on the edges 01, 03 and 04. This system is
  // a rhs system and faces should point outwards the block following the
  // right-hand rule in order to produce a particular block order we need to
  // find out the origin first. Here we follow the convention that our system is
  // (x', y', z) where x' = slope direction and y' = x' x z So the origin vertex
  // will the one with the smallest x' and y'
  auto x_ = hermes::normalize(slope_direction);
  auto y_ = hermes::normalize(x_.left());

  auto localCoordinate = [&](const hermes::point3 &vertex) -> hermes::point2 {
    hermes::vec2 v{vertex.x - slope_point.x, vertex.y - slope_point.y};
    return {
        hermes::dot(v, x_), // x'
        hermes::dot(v, y_)  // y'
    };
  };

  auto findOriginIndex = [&](const std::vector<size_t> &indices) -> size_t {
    size_t origin_index = 0;
    auto origin_point = localCoordinate(input_mesh.vertices[indices[0]]);
    for (size_t i = 1; i < indices.size(); ++i) {
      auto p = localCoordinate(input_mesh.vertices[indices[i]]);
      if (p.x < origin_point.x ||
          (p.x == origin_point.x && p.y < origin_point.y)) {
        origin_index = i;
        origin_point = p;
      }
    }
    return origin_index;
  };

  // 5 - re-order face vertices to make the first the origin
  // //////////////////////////////////////////////////////// we just need to
  // iterate over the bottom surface
  for (size_t face_id = 0; face_id < faces.size(); ++face_id) {
    // find origin vertex index
    size_t origin_index = findOriginIndex(input_mesh.faces[face_id]);
    // shift order to put origin in the first
    std::vector<size_t> new_face_order;
    for (size_t i = 0; i < input_mesh.faces[face_id].size(); ++i)
      new_face_order.emplace_back(
          input_mesh.faces[face_id][(origin_index + i) %
                                    input_mesh.faces[face_id].size()]);
    input_mesh.faces[face_id] = new_face_order;
  }

  // there are two modes for combining blocks: face matching and face merging

  // Face Matching
  // //////////////////////////////////////////////////////////////////////////////////////////////////
  // face matching is less flexible, as vertices must match
  // iterate over original faces (as the extrusion keeps them on the beginning
  // of the mesh)
  if (config.face_matching) {
    for (size_t face_id = 0; face_id < faces.size(); ++face_id) {
      // for each pair (original face, extrusion face), generate a block
      // we know that our faces are properly ordered
      HERMES_ASSERT(face_id < input_mesh.faces.size());
      const auto &face_vertices = input_mesh.faces[face_id];
      BlockMeshDescription::Block block;
      // block vertices
      block.hex = {
          face_vertices[0], // 0
          face_vertices[3], // 1
          face_vertices[2], // 2
          face_vertices[1], // 3

          face_vertices[0] + vertices.size(), // 4
          face_vertices[3] + vertices.size(), // 5
          face_vertices[2] + vertices.size(), // 6
          face_vertices[1] + vertices.size(), // 7
      };
      block.resolution = config.resolution;
      block.grading = config.grading;
      block_mesh.pushBlock(block);

      // boundaries
      // from construction, we know that the top and bottom faces are boundary
      // faces and their indices are 4 and 5
      block_mesh.addFaceToBoundary(bottom_patch_name, block.getFace(4));
      block_mesh.addFaceToBoundary(top_patch_name, block.getFace(5));
      // for the rest of faces, we need to make sure this block is a boundary
      // block
      for (size_t i = 0; i < 4; ++i) {
        auto face = block.getFace(i);
        // first we need to detect if this face is boundary
        // the way to do this is to get which edge-vertices it contains and
        // check if it appears in the edge-face map
        for (size_t j = 0; j < 4; ++j) {
          // since the order is not guaranteed, we need to find if any sequence
          // of vertices appears in the surface face by retrieving its key
          auto edge_key = edgeKey(face[j], face[(j + 1) % 4]);
          auto it = edge_map.find(edge_key);
          if (it != edge_map.end() && it->second.second < 0) {
            // any is sufficient (because we know that the edge_map was
            // constructed from the original faces only)
            auto normal = faceNormal(input_mesh.vertices, face);
            auto label = labelDirection(normal, wall_patch_names);
            if (label.empty()) {
              bad_faces.faces.emplace_back(face);
              HERMES_LOG_ERROR("no label for face with normal {}", normal);
            } else
              block_mesh.addFaceToBoundary(label, face);
          }
        }
      }
    }
  } else {
    HERMES_NOT_IMPLEMENTED;
  }

  saveOBJ(input_mesh, "/home/filipecn/Desktop/extrude.obj");
  saveOBJ(bad_faces, "/home/filipecn/Desktop/bad_faces.obj");

  block_mesh.info();

  return std::move(block_mesh);
}

} // namespace psa_anim
