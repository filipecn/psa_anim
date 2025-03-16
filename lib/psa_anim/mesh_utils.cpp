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
///\file mesh_utils.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-11-30
///
///\brief

#include <hermes/common/index.h>
#include <hermes/data_structures/cartesian_hash.h>
#include <hermes/geometry/bbox.h>
#include <hermes/geometry/queries.h>
#include <hermes/geometry/vector.h>
#include <hermes/logging/logging.h>
#include <hermes/numeric/numeric.h>
#include <mesh_utils.h>
#include <queue>
#include <set>
#define TINYOBJLOADER_IMPLEMENTATION
#include <microstl.h>
#include <tiny_obj_loader.h>

namespace psa_anim {

std::unordered_map<size_t, std::pair<i64, i64>>
computeFaceTopology(const std::vector<std::vector<size_t>> &faces,
                    const std::vector<std::vector<size_t>> &cells,
                    const std::set<size_t> &region) {
  std::unordered_map<size_t, std::pair<i64, i64>> face_cells_map;
  if (region.empty()) {
    for (size_t cell_id = 0; cell_id < cells.size(); ++cell_id)
      for (const auto &face_id : cells[cell_id]) {
        auto it = face_cells_map.find(face_id);
        if (it == face_cells_map.end())
          face_cells_map[face_id] = std::make_pair(cell_id, -1);
        else if (it->second.second >= 0) {
          HERMES_LOG_ERROR("face shared by more than 2 cells.");
          return {};
        } else {
          it->second.second = static_cast<i64>(cell_id);
          if (it->second.first < it->second.second)
            std::swap(it->second.first, it->second.second);
        }
      }
  } else {

    for (size_t cell_id : region)
      for (const auto &face_id : cells[cell_id]) {
        auto it = face_cells_map.find(face_id);
        if (it == face_cells_map.end())
          face_cells_map[face_id] = std::make_pair(cell_id, -1);
        else if (it->second.second >= 0) {
          HERMES_LOG_ERROR("face shared by more than 2 cells.");
          return {};
        } else {
          it->second.second = static_cast<i64>(cell_id);
          if (it->second.first < it->second.second)
            std::swap(it->second.first, it->second.second);
        }
      }
  }
  return std::move(face_cells_map);
}

std::vector<size_t> getNeighbours(
    const std::unordered_map<size_t, std::pair<i64, i64>> &face_topology,
    const std::vector<std::vector<size_t>> &cells, size_t cell_id,
    size_t ring_size) {
  HERMES_ASSERT(cell_id < cells.size());

  std::unordered_map<size_t, size_t> visited;
  std::queue<size_t> q;
  q.push(cell_id);
  visited[cell_id] = 0;
  while (!q.empty()) {
    auto current_cell_id = q.front();
    q.pop();
    HERMES_ASSERT(visited.count(current_cell_id));
    size_t current_level = visited[current_cell_id];
    // get neighbours
    HERMES_ASSERT(current_cell_id < cells.size());
    std::vector<size_t> current_neighbours;
    for (size_t face_id : cells[current_cell_id]) {
      auto it = face_topology.find(face_id);
      HERMES_ASSERT(it != face_topology.end());
      HERMES_ASSERT(it->second.first >= 0);
      if (it->second.first != current_cell_id)
        current_neighbours.emplace_back(it->second.first);
      if (it->second.second != current_cell_id && it->second.second >= 0)
        current_neighbours.emplace_back(it->second.second);
    }
    // filter neighbours
    for (auto nid : current_neighbours) {
      if (visited.count(nid))
        continue;
      if (current_level + 1 > ring_size)
        continue;
      visited[nid] = current_level + 1;
      q.push(nid);
    }
  }
  // neighbours will be at visited
  std::vector<size_t> neighbours;
  for (const auto &v : visited)
    if (v.first != cell_id)
      neighbours.emplace_back(v.first);

  return neighbours;

  for (size_t face_id : cells[cell_id]) {
    auto it = face_topology.find(face_id);
    HERMES_ASSERT(it != face_topology.end())
    HERMES_ASSERT(it->second.first >= 0)
    if (it->second.first != face_id)
      neighbours.emplace_back(it->second.first);
    if (it->second.second != face_id && it->second.second >= 0)
      neighbours.emplace_back(it->second.second);
  }
  return std::move(neighbours);
}

std::vector<size_t>
intersect(const std::unordered_map<size_t, std::pair<i64, i64>> &face_topology,
          const std::vector<std::vector<size_t>> &cells,
          const std::vector<hermes::point3> &centers, size_t cell_id,
          const hermes::bbox3 &region) {
  if (!region.contains(centers[cell_id])) {
    HERMES_LOG_VARIABLES(region, centers[cell_id]);
    HERMES_LOG_ERROR("{} {} {}", cell_id, centers[cell_id], region);
  }
  HERMES_ASSERT(region.contains(centers[cell_id]));
  std::set<size_t> ans;
  std::queue<size_t> q;
  q.push(cell_id);
  ans.insert(cell_id);
  while (!q.empty()) {
    auto current_cell_id = q.front();
    q.pop();
    HERMES_ASSERT(ans.count(current_cell_id));
    // get neighbours
    HERMES_ASSERT(current_cell_id < cells.size());
    std::vector<size_t> current_neighbours;
    for (size_t face_id : cells[current_cell_id]) {
      auto it = face_topology.find(face_id);
      HERMES_ASSERT(it != face_topology.end());
      HERMES_ASSERT(it->second.first >= 0);
      if (it->second.first != current_cell_id)
        current_neighbours.emplace_back(it->second.first);
      if (it->second.second != current_cell_id && it->second.second >= 0)
        current_neighbours.emplace_back(it->second.second);
    }
    // filter neighbours
    for (auto nid : current_neighbours) {
      if (ans.count(nid) || !region.contains(centers[nid]))
        continue;
      ans.insert(nid);
      q.push(nid);
    }
  }
  // neighbours will be at visited
  std::vector<size_t> neighbours;
  for (auto a : ans)
    neighbours.emplace_back(a);

  return neighbours;
}

std::vector<size_t>
getKnn(const std::unordered_map<size_t, std::pair<i64, i64>> &face_topology,
       const std::vector<std::vector<size_t>> &cells,
       const std::vector<hermes::point3> &centers, size_t cell_id, size_t n) {
  HERMES_ASSERT(cell_id < cells.size());

  auto center = centers[cell_id];
  auto cmp = [&](size_t a, size_t b) {
    return hermes::distance2(centers[a], center) >
           hermes::distance2(centers[b], center);
  };
  std::priority_queue<size_t, std::vector<size_t>, decltype(cmp)> pq(cmp);

  std::unordered_map<size_t, size_t> visited;
  std::queue<size_t> q;
  q.push(cell_id);
  visited[cell_id] = 0;
  size_t last_pushed_level = 0;
  while (!q.empty()) {
    auto current_cell_id = q.front();
    q.pop();
    pq.push(current_cell_id);
    if (pq.size() > n * 1.5)
      break;
    HERMES_ASSERT(visited.count(current_cell_id));
    size_t current_level = visited[current_cell_id];
    // get neighbours
    HERMES_ASSERT(current_cell_id < cells.size());
    std::vector<size_t> current_neighbours;
    for (size_t face_id : cells[current_cell_id]) {
      auto it = face_topology.find(face_id);
      HERMES_ASSERT(it != face_topology.end());
      HERMES_ASSERT(it->second.first >= 0);
      if (it->second.first != current_cell_id)
        current_neighbours.emplace_back(it->second.first);
      if (it->second.second != current_cell_id && it->second.second >= 0)
        current_neighbours.emplace_back(it->second.second);
    }
    // filter neighbours
    for (auto nid : current_neighbours) {
      if (visited.count(nid))
        continue;
      visited[nid] = current_level + 1;
      q.push(nid);
    }
  }
  // neighbours will be at visited
  std::vector<size_t> neighbours;
  while (!pq.empty() && neighbours.size() < n) {
    neighbours.emplace_back(pq.top());
    pq.pop();
  }

  return neighbours;
}

real_t directionAlignment(const hermes::vec3 &a, const hermes::vec3 &b) {
  return hermes::dot(hermes::normalize(a), hermes::normalize(b));
}

bool directionIsTheSame(const hermes::vec3 &a, const hermes::vec3 &b) {
  return directionAlignment(a, b) > 0;
}

hermes::point3 faceCenter(const std::vector<hermes::point3> &vertices,
                          const std::vector<size_t> &face_vertices) {
  hermes::point3 numerator;
  for (auto vertex : face_vertices) {
    HERMES_ASSERT(vertex < vertices.size());
    numerator += (hermes::vec3)vertices[vertex];
  }
  HERMES_ASSERT(!face_vertices.empty());
  return numerator / (f32)face_vertices.size();
}

hermes::point3 cellCenter(const std::vector<hermes::point3> &vertices,
                          const std::vector<std::vector<size_t>> &faces,
                          const std::vector<size_t> &cell_faces) {
  std::set<size_t> cell_vertices;
  for (size_t cell_face : cell_faces) {
    HERMES_ASSERT(cell_face < faces.size());
    for (size_t vertex_id : faces[cell_face])
      cell_vertices.insert(vertex_id);
  }
  hermes::point3 numerator;
  for (auto vertex : cell_vertices)
    numerator += (hermes::vec3)vertices[vertex];
  HERMES_ASSERT(!cell_vertices.empty());
  return numerator / (f32)cell_vertices.size();
}

std::vector<hermes::point3>
cellCenters(const std::vector<hermes::point3> &vertices,
            const std::vector<std::vector<size_t>> &faces,
            const std::vector<std::vector<size_t>> &cells) {
  std::vector<hermes::point3> centers;
  for (const auto &cell_faces : cells)
    centers.emplace_back(cellCenter(vertices, faces, cell_faces));
  return std::move(centers);
}

hermes::vec3 faceNormal(const std::vector<hermes::point3> &vertices,
                        const std::vector<size_t> &face_vertices) {
  HERMES_ASSERT(face_vertices.size() > 2);
  // TODO naive (using just the first triangle
  return hermes::normalize(
      hermes::cross(vertices[face_vertices[1]] - vertices[face_vertices[0]],
                    vertices[face_vertices[2]] - vertices[face_vertices[1]]));
}

void fixBoundaryFaces(const std::vector<hermes::point3> &vertices,
                      const std::vector<std::vector<size_t>> &cells,
                      std::vector<std::vector<size_t>> &faces) {
  // compute data
  auto face_cells_map = computeFaceTopology(faces, cells);
  auto cell_centers = cellCenters(vertices, faces, cells);
  HERMES_ASSERT(!face_cells_map.empty());
  // flip faces pointing inwards the mesh
  for (const auto &face : face_cells_map) {
    HERMES_ASSERT(face.second.first >= 0);
    if (face.second.second < 0) {
      // TODO naive (testing normal with cell center)
      // normal should point away from center
      auto face_center = faceCenter(vertices, faces[face.first]);
      auto center_direction = cell_centers[face.second.first] - face_center;
      auto normal = faceNormal(vertices, faces[face.first]);
      if (hermes::dot(center_direction, normal) > 0) {
        // normal and center pointing on the same direction, face must be
        // flipped
        auto face_copy = faces[face.first];
        for (size_t i = 0; i < face_copy.size(); ++i)
          face_copy[i] = faces[face.first][face_copy.size() - i - 1];
      }
    }
  }
}

std::vector<std::pair<size_t, size_t>> markBoundary(
    const std::vector<hermes::point3> &vertices,
    std::vector<std::vector<size_t>> &cells,
    std::vector<std::vector<size_t>> &faces,
    const std::vector<std::vector<hermes::vec3>> &boundary_directions) {
  auto normalized_boundary_directions = boundary_directions;
  for (auto &vb : normalized_boundary_directions)
    for (auto &v : vb)
      v.normalize();
  std::unordered_map<size_t, size_t> face_boundary_map;
  // detect boundary faces and separate them into their respective set
  auto face_cells_map = computeFaceTopology(faces, cells);
  size_t boundary_face_count = 0;
  for (size_t face_id = 0; face_id < faces.size(); ++face_id) {
    auto it = face_cells_map.find(face_id);
    HERMES_ASSERT(it != face_cells_map.end());
    if (it->second.second < 0) {
      boundary_face_count++;
      auto face_normal = faceNormal(vertices, faces[face_id]);
      // find out which boundary this face belongs
      auto largest_dot_value = 0.f;
      size_t largest_id = 0;
      for (size_t boundary_id = 0;
           boundary_id < normalized_boundary_directions.size(); boundary_id++) {
        for (const auto &direction :
             normalized_boundary_directions[boundary_id]) {
          auto dot_value = hermes::dot(face_normal, direction);
          if (dot_value > largest_dot_value) {
            largest_dot_value = dot_value;
            largest_id = boundary_id;
          }
        }
      }
      // label face
      face_boundary_map[face_id] = largest_id;
    }
  }
  HERMES_LOG_VARIABLE(boundary_face_count);
  // sort faces
  std::vector<size_t> face_indices(faces.size());
  for (size_t i = 0; i < face_indices.size(); ++i)
    face_indices[i] = i;
  std::sort(face_indices.begin(), face_indices.end(), [&](size_t a, size_t b) {
    auto a_it = face_boundary_map.find(a);
    auto b_it = face_boundary_map.find(b);
    if (a_it == face_boundary_map.end() && b_it != face_boundary_map.end())
      return true;
    if (b_it == face_boundary_map.end() && a_it != face_boundary_map.end())
      return false;
    if (a_it == face_boundary_map.end() && b_it == face_boundary_map.end())
      return a < b;
    return a_it->second < b_it->second;
  });
  // retrieve new indices for faces
  std::unordered_map<size_t, size_t> face_index_map;
  for (size_t i = 0; i < face_indices.size(); ++i)
    face_index_map[face_indices[i]] = i;
  // re-label cell face indices
  for (auto &cell_faces : cells)
    for (auto &face_id : cell_faces) {
      HERMES_ASSERT(face_index_map.count(face_id));
      face_id = face_index_map[face_id];
    }
  // re-order faces
  // TODO re-order inline instead of copying the entire face array
  // also get boundary patches
  std::vector<std::pair<size_t, size_t>> patches;
  for (size_t i = 0; i < boundary_directions.size(); ++i)
    patches.emplace_back(faces.size(), 0);
  //
  std::vector<std::vector<size_t>> new_faces(faces.size());
  for (size_t face_id = 0; face_id < face_indices.size(); face_id++) {
    HERMES_ASSERT(face_index_map.count(face_id));
    new_faces[face_index_map[face_id]] = faces[face_id];
    auto it = face_boundary_map.find(face_id);
    if (it != face_boundary_map.end()) {
      patches[it->second].first =
          std::min(patches[it->second].first, face_index_map[face_id]);
      patches[it->second].second++;
    }
  }
  faces = new_faces;
  return patches;
}

std::map<std::string, FaceMesh> splitBoundary(
    const FaceMesh &mesh,
    const std::map<std::string, std::vector<hermes::vec3>> &boundary_directions,
    const std::map<std::string, std::vector<std::pair<size_t, size_t>>>
        &boundary_ranges) {
  auto normalized_boundary_directions = boundary_directions;
  for (auto &vb : normalized_boundary_directions)
    for (auto &v : vb.second)
      v.normalize();

  HERMES_LOG("++++++++++++++++++++++++++++++++++++++++++");
  HERMES_LOG("Splitting mesh into boundary patches:");
  for (const auto &boundary : boundary_directions) {
    HERMES_LOG("boundary patch: {}", boundary.first);
    HERMES_LOG_ARRAY(boundary.second);
  }

  // label faces
  std::unordered_map<size_t, std::string> face_boundary_map;
  for (size_t face_id = 0; face_id < mesh.faces.size(); ++face_id) {
    std::string largest_id;
    // check against index ranges
    for (const auto &boundary_range : boundary_ranges)
      for (const auto r : boundary_range.second)
        if (face_id >= r.first && face_id < r.second)
          largest_id = boundary_range.first;
    if (largest_id.empty()) {
      auto face_normal = faceNormal(mesh.vertices, mesh.faces[face_id]);
      // find out which boundary this face belongs
      auto largest_dot_value = 0.f;
      for (const auto &boundary : normalized_boundary_directions) {
        for (const auto &direction : boundary.second) {
          auto dot_value = hermes::dot(face_normal, direction);
          if (dot_value > largest_dot_value) {
            largest_dot_value = dot_value;
            largest_id = boundary.first;
          }
        }
      }
    }
    // label face
    HERMES_ASSERT(!largest_id.empty());
    face_boundary_map[face_id] = largest_id;
  }

  // split faces into patches (renaming the vertices)
  std::map<std::string, FaceMesh> patches;
  // old vertex id -> new vertex id
  std::map<std::string, std::unordered_map<size_t, size_t>> patch_vertex_id;

  // split
  for (size_t face_id = 0; face_id < mesh.faces.size(); ++face_id) {
    // get patch for this face
    HERMES_ASSERT(face_id < face_boundary_map.size());
    auto patch_name = face_boundary_map[face_id];
    std::vector<size_t> face;
    for (auto old_vertex_id : mesh.faces[face_id]) {
      if (!patch_vertex_id[patch_name].count(old_vertex_id)) {
        patch_vertex_id[patch_name][old_vertex_id] =
            patches[patch_name].vertices.size();
        patches[patch_name].vertices.emplace_back(mesh.vertices[old_vertex_id]);
      }
      face.emplace_back(patch_vertex_id[patch_name][old_vertex_id]);
    }
    patches[patch_name].faces.emplace_back(face);
  }

  return patches;
}

std::pair<size_t, size_t> edgeKey(size_t va, size_t vb) {
  return {std::min(va, vb), std::max(va, vb)};
}

std::map<std::pair<size_t, size_t>, std::pair<i64, i64>>
computeEdgeTopology(const std::vector<std::vector<size_t>> &faces) {
  std::map<std::pair<size_t, size_t>, std::pair<i64, i64>> edge_face_map;
  // for each facet, register the new edges
  for (size_t t = 0; t < faces.size(); ++t) {
    const auto &face_vertices = faces[t];
    // iterate edges
    for (int i = 0; i < face_vertices.size(); ++i) {
      auto key = edgeKey(face_vertices[i],
                         face_vertices[(i + 1) % face_vertices.size()]);
      if (!edge_face_map.count(key))
        edge_face_map[key] = {t, -1};
      else {
        HERMES_ASSERT(edge_face_map[key].second == -1);
        edge_face_map[key].second = t;
      }
    }
  }
  return std::move(edge_face_map);
}

std::vector<size_t> getSurfaceNeighbours(
    const std::map<std::pair<size_t, size_t>, std::pair<i64, i64>>
        &edge_topology,
    const std::vector<std::vector<size_t>> &faces, size_t face_id) {
  HERMES_ASSERT(face_id < faces.size())
  std::vector<size_t> neighbours;
  const auto &face = faces[face_id];
  for (size_t i = 0; i < face.size(); ++i) {
    auto edge_key = edgeKey(face[i], face[(i + 1) % face.size()]);
    auto it = edge_topology.find(edge_key);
    HERMES_ASSERT(it != edge_topology.end())
    HERMES_ASSERT(it->second.first >= 0)
    if (it->second.first != face_id)
      neighbours.emplace_back(it->second.first);
    if (it->second.second != face_id && it->second.second >= 0)
      neighbours.emplace_back(it->second.second);
  }
  return neighbours;
}

void flipOrder(std::vector<size_t> &v) {
  for (size_t i = 0; i < v.size() / 2; ++i)
    std::swap(v[i], v[v.size() - 1 - i]);
}

FaceMesh genXYAlignedMesh(const FaceMesh &input, const hermes::point2 &origin,
                          float cell_size) {
  // setup grid dimensions
  hermes::bbox2 xy_bounds;
  for (const auto &v : input.vertices)
    xy_bounds = hermes::make_union(xy_bounds, v.xy());
  hermes::vec2 size(xy_bounds.size(0), xy_bounds.size(1));
  HERMES_LOG_VARIABLE(xy_bounds);
  HERMES_LOG_VARIABLE(size);

  auto quadMesh = extractGrid(input, {1, 0}, origin, size, cell_size, false);
  FaceMesh m;
  m.vertices = quadMesh.vertices;
  for (const auto &face : quadMesh.faces) {
    m.faces.push_back({face[0], face[1], face[2]});
    m.faces.push_back({face[0], face[2], face[3]});
  }
  return m;
}

FaceMesh extrudeZ(const FaceMesh &surface, real_t height,
                  hermes::vec2 slope_direction, hermes::point2 slope_point,
                  bool extrusion_surface_is_inclined,
                  bool generate_boundary_faces,
                  bool boundary_faces_are_triangles) {

  // prepare input
  slope_direction = hermes::normalize(slope_direction);

  HERMES_LOG("Extruding surface parameters");
  HERMES_LOG("==================================================");
  HERMES_LOG_VARIABLE(height);
  HERMES_LOG_VARIABLE(slope_direction);
  HERMES_LOG_VARIABLE(slope_point);
  HERMES_LOG_VARIABLE((int)extrusion_surface_is_inclined);
  HERMES_LOG_VARIABLE((int)generate_boundary_faces);
  HERMES_LOG_VARIABLE((int)boundary_faces_are_triangles);
  HERMES_LOG("==================================================");
  // data
  FaceMesh mesh;
  auto edge_face_map = computeEdgeTopology(surface.faces);

  // The extrusion algorithm follows the steps:
  // 0 - reserve new memory for the new mesh and copy the surface
  // 1 - fix terrain and top normals
  // 2 - compute surface bounds and define the top height
  // 3 - compute the top plane
  // 4 - project the duplicates into the top plane
  // 5 - generate boundary faces

  // 0 reserve memory
  // ///////////////////////////////////////////////////////////////////////////////////////////////
  HERMES_LOG("reserving memory for vertices:{} faces:{}",
             surface.vertices.size() * 2, surface.faces.size() * 2);
  // copy the original surface and duplicate the vertices
  mesh.vertices.resize(surface.vertices.size() * 2);
  // we also duplicate faces (and compute the new indices)
  mesh.faces.resize(surface.faces.size() * 2);
  for (size_t i = 0; i < surface.vertices.size(); ++i)
    mesh.vertices[i] = mesh.vertices[i + surface.vertices.size()] =
        surface.vertices[i];
  for (size_t i = 0; i < surface.faces.size(); ++i) {
    mesh.faces[i] = surface.faces[i];
    mesh.faces[i + surface.faces.size()] = surface.faces[i];
    for (auto &v : mesh.faces[i + surface.faces.size()])
      v += surface.vertices.size();
    // fix face normals
    if (!directionIsTheSame(hermes::vec3(0, 0, -1),
                            faceNormal(mesh.vertices, mesh.faces[i])))
      flipOrder(mesh.faces[i]);
    else
      flipOrder(mesh.faces[i + surface.faces.size()]);
  }

  // 1-3 top patch
  // //////////////////////////////////////////////////////////////////////////////////////////////////
  // the bounds must be aligned to the slope direction
  // we project all points into the slope line. The projection distance computes
  // the projection of a vertex onto the axis line and the distance of the
  // projected point to the line start
  auto projectionDistance = [&](const hermes::point3 &vertex) -> real_t {
    hermes::vec2 v{vertex.x - slope_point.x, vertex.y - slope_point.y};
    return hermes::dot(v, slope_direction);
  };

  // read left and right here in the sense of the slope direction, which points
  // to the "right"
  real_t left_most_projection = projectionDistance(surface.vertices[0]);
  real_t right_most_projection = left_most_projection;
  // in order to find the top plane, we need to get its descent by retrieving
  // the bounding z values
  real_t min_z = surface.vertices[0].z;
  real_t max_z = surface.vertices[0].z;

  for (const auto &vertex : surface.vertices) {
    auto current_projection = projectionDistance(vertex);
    left_most_projection = std::min(left_most_projection, current_projection);
    right_most_projection = std::max(right_most_projection, current_projection);
    min_z = std::min(min_z, vertex.z);
    max_z = std::max(max_z, vertex.z);
  }

  // from now on lets call the new x' axis as the slope_direction, and the new
  // y' axis as the cross(x', z) assuming that the lower height is at the
  // right-most point the plane has no rotation in x', only in y', and can be
  // computed from the heights of the two extremes:
  real_t middle_axis_distance =
      (left_most_projection + right_most_projection) * 0.5f;
  hermes::point2 m_axis_a =
      slope_point + left_most_projection * slope_direction;
  hermes::point2 m_axis_b =
      slope_point + right_most_projection * slope_direction;
  hermes::vec2 y_axis = slope_direction.left();
  // now we include z information into the slope direction
  hermes::point3 a(m_axis_a, max_z + height * 0.5);
  hermes::point3 b(m_axis_b, min_z + height);
  hermes::vec3 descenting_slope_direction = hermes::normalize(b - a);
  if (!extrusion_surface_is_inclined)
    // then top is aligned to XY
    descenting_slope_direction = hermes::vec3(descenting_slope_direction.x,
                                              descenting_slope_direction.y, 0)
                                     .normalized();

  HERMES_LOG_VARIABLE(min_z);
  HERMES_LOG_VARIABLE(max_z);
  HERMES_LOG_VARIABLE(a);
  HERMES_LOG_VARIABLE(b);
  HERMES_LOG_VARIABLE(left_most_projection);
  HERMES_LOG_VARIABLE(right_most_projection);
  HERMES_LOG_VARIABLE(m_axis_a);
  HERMES_LOG_VARIABLE(m_axis_b);
  HERMES_LOG_VARIABLE(y_axis);
  HERMES_LOG_VARIABLE(middle_axis_distance);
  HERMES_LOG_VARIABLE(descenting_slope_direction);

  // 4 - project the duplicates into the top plane
  // ////////////////////////////////////////////////////////////////// now we
  // project the points into the top plane (in z direction)
  for (size_t i = surface.vertices.size(); i < mesh.vertices.size(); ++i) {
    auto &vertex = mesh.vertices[i];
    auto current_projection = projectionDistance(vertex);
    // the projection goes like this:
    vertex.z = (a + current_projection * descenting_slope_direction).z;
  }

  if (!generate_boundary_faces)
    return std::move(mesh);

  // 6 - generate boundary faces
  // ////////////////////////////////////////////////////////////////////////////////////
  // at this point, we have the bottom and top surfaces
  // and we also have the edge map for the bottom surface
  // the edges were duplicated and the vertices too
  for (const auto &edge : edge_face_map) {
    // a boundary edge does not has the second neighbour
    HERMES_ASSERT(edge.second.first >= 0);
    // skip internal faces
    if (edge.second.second >= 0)
      continue;
    // find out which patch this edge belongs to
    auto edge_a = edge.first.first;
    auto edge_b = edge.first.second;
    // the edge_a -> edge_b order comes from the key, not the original order
    // we must flip if not in the same order of the bottom surface
    if (!mesh.edgeIsSorted(edge.second.first, edge_a, edge_b))
      std::swap(edge_a, edge_b);
    // now we know that for the same face in the top patch the edge is flipped:
    // edge_b -> edge_a in order to create a normal pointing outside the mesh,
    // we need to follow the edges order
    if (boundary_faces_are_triangles) {
      mesh.faces.push_back({edge_a, edge_b, edge_a + surface.vertices.size()});
      mesh.faces.push_back({edge_b, edge_b + surface.vertices.size(),
                            edge_a + surface.vertices.size()});
    } else {
      mesh.faces.push_back({edge_a, edge_b, edge_b + surface.vertices.size(),
                            edge_a + surface.vertices.size()});
    }
  }

  // return
  return std::move(mesh);
}

hermes::Result<FaceMesh> loadOBJ(const hermes::Path &path) {
  // return
  FaceMesh mesh;

  tinyobj::ObjReaderConfig reader_config;
  reader_config.triangulate = false;
  tinyobj::ObjReader reader;

  if (!reader.ParseFromFile(path.fullName(), reader_config)) {
    if (!reader.Error().empty())
      HERMES_LOG_ERROR(reader.Error().c_str());
    return hermes::Result<FaceMesh>::error(HeResult::INVALID_INPUT);
  }

  if (!reader.Warning().empty())
    HERMES_LOG_WARNING(reader.Warning().c_str());

  auto &attrib = reader.GetAttrib();
  auto &shapes = reader.GetShapes();

  // retrieve vertices
  auto n_vertices = attrib.vertices.size() / 3;
  for (size_t i = 0; i < n_vertices; ++i)
    mesh.vertices.emplace_back(attrib.vertices[i * 3 + 0],
                               attrib.vertices[i * 3 + 1],
                               attrib.vertices[i * 3 + 2]);

  // Loop over shapes
  for (const auto &shape : shapes) {
    // Loop over faces(polygon)
    size_t index_offset = 0;
    for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) {
      auto fv = size_t(shape.mesh.num_face_vertices[f]);
      std::vector<size_t> face_vertices;
      // Loop over vertices in the face.
      for (size_t v = 0; v < fv; v++) {
        // access to vertex
        tinyobj::index_t idx = shape.mesh.indices[index_offset + v];
        face_vertices.emplace_back(idx.vertex_index);
      }
      index_offset += fv;
      mesh.faces.emplace_back(face_vertices);
    }
    break;
  }

  return hermes::Result<FaceMesh>(mesh);
}

bool saveOBJ(const FaceMesh &mesh, const hermes::Path &path, bool wireframe) {
  hermes::Str s;
  for (const auto &v : mesh.vertices)
    s.appendLine("v ", v.x, " ", v.y, " ", v.z);
  for (const auto &f : mesh.faces) {
    if (wireframe)
      s.append("l");
    else
      s.append("f");
    for (auto ff : f)
      s.append(" ", ff + 1);
    s.appendLine("");
  }
  hermes::FileSystem::writeFile(path, s.str());
  return true;
}

bool saveSTL(const std::map<std::string, FaceMesh> &mesh,
             const hermes::Path &path) {
  HERMES_LOG("writing into {}.", path.fullName());
  hermes::Str str;
  bool split_files = true;
  if (path.isFile() || path.hasExtension())
    split_files = false;
  auto write = [&](const std::string &name) {
    auto it = mesh.find(name);
    HERMES_ASSERT(it != mesh.end());
    // convert to microstl
    microstl::FVMesh micro_mesh;
    for (const auto &v : it->second.vertices)
      micro_mesh.vertices.push_back({v.x, v.y, v.z});
    for (const auto &f : it->second.faces) {
      HERMES_ASSERT(f.size() == 3);
      auto n = faceNormal(it->second.vertices, f);
      micro_mesh.facets.push_back(
          {.v1 = f[0], .v2 = f[1], .v3 = f[2], .n = {n.x, n.y, n.z}});
    }
    microstl::FVMeshProvider provider(micro_mesh);
    provider.ascii = true;
    std::string s;
    microstl::Writer::writeStlBuffer(s, provider);
    if (split_files) {
      auto filename = path / (name + ".stl");
      hermes::FileSystem::writeFile(
          filename, hermes::Str::regex::replace(s, provider.getName(), name));
    } else
      str.append(hermes::Str::regex::replace(s, provider.getName(), name));
  };

  for (const auto &patch : mesh)
    write(patch.first);

  if (!split_files)
    hermes::FileSystem::writeFile(path, str.c_str());

  return true;
}

hermes::Result<FaceMesh> parseSTL(const std::string &buffer) {
  HERMES_LOG("calling microstl...");
  microstl::MeshReaderHandler mesh_handler;
  auto result = microstl::Reader::readStlBuffer(buffer.c_str(), buffer.size(),
                                                mesh_handler);
  if (result != microstl::Result::Success) {
    HERMES_LOG_ERROR("stl load error: {}", microstl::getResultString(result));
    return hermes::Result<FaceMesh>::error(HeResult::INVALID_INPUT);
  }
  auto stl_surface = microstl::deduplicateVertices(mesh_handler.mesh);
  HERMES_LOG_VARIABLE(stl_surface.vertices.size());
  HERMES_LOG_VARIABLE(stl_surface.facets.size());
  FaceMesh mesh;
  for (auto &v : stl_surface.vertices)
    mesh.vertices.emplace_back(v.x, v.y, v.z);
  for (auto &f : stl_surface.facets)
    mesh.faces.push_back({f.v1, f.v2, f.v3});
  return hermes::Result<FaceMesh>(mesh);
}

hermes::Result<FaceMesh> loadSTL(const hermes::Path &path) {
  auto buffer = hermes::FileSystem::readFile(path);
  return parseSTL(buffer);
}

hermes::Result<FaceMesh> loadSTL(const hermes::Path &path,
                                 const std::string &solid) {
  auto s = hermes::FileSystem::readFile(path);
  // filter patch from stl
  size_t i = 0;
  size_t j = 0;
  size_t n = s.size();
  bool found = false;
  while (i < n && !found) {
    if (s[i] == 's' && s.substr(i, 5) == "solid") {
      // check patch name
      j = i + 5;
      while (s[j] == ' ')
        j++;
      if (s.substr(j, solid.size()) == solid) {
        // found!
        while (j < n) {
          if (s[j] == 'e' && s.substr(j, 8) == "endsolid") {
            j += 9;
            found = true;
            break;
          }
          j++;
        }
      }
    }
    ++i;
  }
  if (!found)
    return hermes::Result<FaceMesh>::error(HeResult::INVALID_INPUT);
  HERMES_ASSERT(i > 0 && j > i);
  return parseSTL(s.substr(i - 1, j - i));
}

FaceMesh convertToTri(const FaceMesh &input) {
  FaceMesh mesh;
  mesh.vertices = input.vertices;

  for (size_t face_id = 0; face_id < input.faces.size(); ++face_id) {
    const auto &face = input.faces[face_id];
    size_t face_size = face.size();
    HERMES_ASSERT(face_size >= 3);
    if (face_size == 3) {
      mesh.faces.push_back(input.faces[face_id]);
    } else {
      auto center_id = mesh.vertices.size();
      mesh.vertices.emplace_back(input.faceCenter(face_id));
      for (size_t v = 0; v < face_size; ++v) {
        mesh.faces.push_back({face[v], face[(v + 1) % face_size], center_id});
      }
    }
  }
  return mesh;
}

FaceMesh convertToQuad(const FaceMesh &input) {
  // return
  FaceMesh mesh;
  // input must be a triangular mesh
  for (const auto &face : input.faces)
    HERMES_ASSERT(face.size() == 3);
  HERMES_ASSERT(input.faces.size() % 2 == 0);
  // the vertices keep the same
  mesh.vertices = input.vertices;
  auto total_faces = input.faces.size() / 2;
  for (size_t i = 0; i < total_faces; ++i) {
    // merge faces i and i + 1
    // both must have the same normal
    HERMES_ASSERT(
        directionIsTheSame(faceNormal(input.vertices, input.faces[i]),
                           faceNormal(input.vertices, input.faces[i + 1])));
    // find out the diagonal edge, and remove it
    std::vector<size_t> face;
    for (size_t j = 0; j < 3; ++j) {
      auto edge = edgeKey(input.faces[i][j], input.faces[i][(j + 1) % 3]);
      // find this edge on the next triangle
      for (size_t k = 0; k < 3; ++k) {
        auto e = edgeKey(input.faces[i][k], input.faces[i][(k + 1) % 3]);
        if (edge == e) {
          // there can be no more than one edge shared!
          HERMES_ASSERT(face.empty());
          // that means that the following vertices are the same:
          HERMES_ASSERT(input.faces[i][j] == input.faces[i][(k + 1) % 3]);
          HERMES_ASSERT(input.faces[i][k] == input.faces[i][(j + 1) % 3]);
          // which also helps us construct the final quad:
          face.emplace_back(input.faces[i][k]);
          face.emplace_back(input.faces[i][(j + 2) % 3]);
          face.emplace_back(input.faces[i][j]);
          face.emplace_back(input.faces[i][(k + 2) % 3]);
        }
      }
    }
    HERMES_ASSERT(!face.empty());
    mesh.faces.emplace_back(face);
  }
  return std::move(mesh);
}

FaceMesh extractGrid(const FaceMesh &input, const hermes::vec2 &axis,
                     const hermes::point2 &origin, const hermes::vec2 &size,
                     real_t cell_size, bool twoDim) {
  // check input mesh type
  for (const auto &face : input.faces) {
    HERMES_ASSERT(face.size() == 3);
    for (auto f : face)
      HERMES_ASSERT(f < input.vertices.size());
  }

  HERMES_LOG("extractGrid:");
  HERMES_LOG_VARIABLE(input.vertices.size());
  HERMES_LOG_VARIABLE(input.faces.size());
  HERMES_LOG_VARIABLE(axis);
  HERMES_LOG_VARIABLE(origin);
  HERMES_LOG_VARIABLE(size);
  HERMES_LOG_VARIABLE(cell_size);
  // prepare search data structure

  // Here we create a grid aligned to an arbitrary axis
  // In its width along the perpendicular direction of the axis, the origin lies
  // in the middle:
  //            |------|------|               ---
  //            |------|------|                |
  //     origin *------*------*  ---> axis     |
  //            |------|------|                |   width
  //            |------|------|               __
  //
  // In other words, a new coordinate system is constructed:
  //   x' = axis
  //   y' = (-axis.y, axis.x)
  auto x_ = hermes::normalize(axis);
  auto y_ = hermes::normalize(x_.left());

  // The algorithm goes as follows
  // 1 - generate grid
  // 2 - compute vertices heights
  // 3 - handle outside faces

  // output mesh
  FaceMesh grid_mesh;
  // inside/outside vertex flag
  std::unordered_map<size_t, bool> vertex_is_inside;

  HERMES_LOG("generating grid...");
  // 1 - generate grid
  // /////////////////////////////////////////////////////////////////////////////////////////////////
  // compute grid resolution
  hermes::size2 res(size.x / cell_size + 1, size.y / cell_size + 1);
  if (twoDim)
    res.height = 2;
  auto flat = [&](const auto &ij) -> size_t { return res.width * ij.j + ij.i; };

  auto half_res = res / 2;
  for (auto ij : hermes::range2(res)) {
    if (ij < res - hermes::index2{1, 1}) {
      grid_mesh.faces.push_back({flat(ij.plus(0, 0)), flat(ij.plus(0, 1)),
                                 flat(ij.plus(1, 1)), flat(ij.plus(1, 0))});
    }
    ij.j -= half_res.height;
    auto point = origin + cell_size * (static_cast<real_t>(ij.i) * x_ +
                                       static_cast<real_t>(ij.j) * y_);
    if (twoDim)
      grid_mesh.vertices.emplace_back(point.x, point.y + cell_size / 2, 0);
    else
      grid_mesh.vertices.emplace_back(point.x, point.y, 0);
  }

  HERMES_LOG("preparing vertices for projection...");
  // 2 - project vertices
  // //////////////////////////////////////////////////////////////////////////////////////////////
  // setup search structure to optimize projection
  std::vector<hermes::point2> terrain_face_positions;
  for (const auto &face : input.faces)
    terrain_face_positions.emplace_back(faceCenter(input.vertices, face).xy());
  auto terrain_face_positions_set =
      hermes::CartesianHashMap2::from(terrain_face_positions, cell_size * 10);

  std::vector<hermes::point2> grid_face_positions;
  for (const auto &face : grid_mesh.faces)
    grid_face_positions.emplace_back(faceCenter(grid_mesh.vertices, face).xy());
  auto grid_face_positions_set =
      hermes::CartesianHashMap2::from(grid_face_positions, cell_size * 10);

  // for each vertex, intersect the input mesh triangles
  HERMES_LOG("projecting vertices...");
  for (size_t vertex_id = 0; vertex_id < grid_mesh.vertices.size();
       ++vertex_id) {
    auto &vertex = grid_mesh.vertices[vertex_id];
    vertex_is_inside[vertex_id] = false;
    auto count = terrain_face_positions_set.search(
        vertex.xy(), cell_size * 4, [&](size_t face_id) {
          HERMES_ASSERT(face_id < input.faces.size());
          const auto &face = input.faces[face_id];
          HERMES_ASSERT(input.vertices.size() > face[0]);
          HERMES_ASSERT(input.vertices.size() > face[1]);
          if (input.vertices.size() <= face[2])
            HERMES_LOG_VARIABLES(face_id, face.size(), input.vertices.size(),
                                 face[2], input.faces[face_id][2]);
          HERMES_ASSERT(input.vertices.size() > face[2]);
          auto a = input.vertices[face[0]];
          auto b = input.vertices[face[1]];
          auto c = input.vertices[face[2]];
          // our ray always goes on the -z direction (terrain)
          hermes::ray3 ray(
              {vertex.x, vertex.y, 10 + std::max(a.z, std::max(b.z, c.z))},
              {0, 0, -1});
          auto t = hermes::GeometricPredicates::intersect(a, b, c, ray);
          if (t) {
            // a hit!
            vertex.z = ray(*t).z;
            vertex_is_inside[vertex_id] = true;
            return false;
          }
          return true;
        });
  }

  HERMES_LOG("extending boundaries...");
  // 3 - handle outside faces

  auto countOutsideVerticesInFace =
      [&](const std::vector<size_t> &face) -> size_t {
    size_t outside_vertex_count = 0;
    for (auto v : face)
      if (!vertex_is_inside[v])
        outside_vertex_count++;
    return outside_vertex_count;
  };

  auto faceIsInside = [&](const std::vector<size_t> &face) -> bool {
    return countOutsideVerticesInFace(face) == 0;
  };

  // here we will propagate the projection outwards the inside domain
  // so we need to get the first layer of outside faces and propagate from there
  // through bfs
  auto grid_mesh_topology = computeEdgeTopology(grid_mesh.faces);
  std::vector<int> visited(grid_mesh.faces.size());
  std::queue<size_t> q;
  for (size_t face_id = 0; face_id < grid_mesh.faces.size(); ++face_id) {
    visited[face_id] = true;
    // filter out inside faces
    if (countOutsideVerticesInFace(grid_mesh.faces[face_id]) == 0)
      continue;
    visited[face_id] = false;
    // if this face neighbours an inside face, add it to the queue
    auto neighbours =
        getSurfaceNeighbours(grid_mesh_topology, grid_mesh.faces, face_id);
    for (auto neighbour_id : neighbours)
      if (faceIsInside(grid_mesh.faces[neighbour_id])) {
        q.push(face_id);
        break;
      }
  }
  // now we propagate
  while (!q.empty()) {
    auto face_id = q.front();
    q.pop();

    if (visited[face_id])
      continue;
    visited[face_id] = true;

    // project vertices based on neighbours
    // get neighbouring inside faces and compute the mean normal
    hermes::point3 plane_center;
    hermes::normal3 plane_normal;
    size_t inside_count = 0;
    auto face_center =
        faceCenter(grid_mesh.vertices, grid_mesh.faces[face_id]).xy();

    size_t count = grid_face_positions_set.search(
        face_center, cell_size * 4, [&](size_t neighbour_face_id) {
          if (neighbour_face_id != face_id && visited[neighbour_face_id]) {
            inside_count++;
            plane_center += (hermes::vec3)faceCenter(
                grid_mesh.vertices, grid_mesh.faces[neighbour_face_id]);
            auto face_normal = faceNormal(grid_mesh.vertices,
                                          grid_mesh.faces[neighbour_face_id]);
            HERMES_ASSERT(hermes::dot(face_normal, hermes::vec3(0, 0, -1)) > 0);
            plane_normal += -faceNormal(grid_mesh.vertices,
                                        grid_mesh.faces[neighbour_face_id]);
          }
          return true;
        });
    HERMES_ASSERT(inside_count)
    plane_center /= static_cast<real_t>(inside_count);
    plane_normal /= static_cast<real_t>(inside_count);
    plane_normal = hermes::normalize(plane_normal);
    hermes::Plane plane(plane_normal, plane_center);
    // for each outside vertex, project it onto the plane
    for (auto v : grid_mesh.faces[face_id])
      if (!vertex_is_inside[v]) {
        hermes::Ray3 ray(grid_mesh.vertices[v] - hermes::vec3(0, 0, 1000),
                         {0, 0, 1});
        auto t = hermes::GeometricPredicates::intersect(plane, ray);
        HERMES_ASSERT(t)
        grid_mesh.vertices[v].z = ray(*t).z;
        vertex_is_inside[v] = true;
      }

    // queue neighbours
    auto neighbours =
        getSurfaceNeighbours(grid_mesh_topology, grid_mesh.faces, face_id);

    for (auto neighbour_id : neighbours)
      if (!visited[neighbour_id])
        q.push(neighbour_id);
  }
  return grid_mesh;
}

void info(const FaceMesh &surface) {
  std::unordered_map<size_t, size_t> face_size;
  for (const auto &face : surface.faces)
    face_size[face.size()]++;
  HERMES_LOG("Surface Info:");
  HERMES_LOG("faces[{}] vertices[{}]", surface.faces.size(),
             surface.vertices.size());
  for (auto &it : face_size)
    HERMES_LOG("#face[{}] = {}", it.first, it.second);
}

} // namespace psa_anim
