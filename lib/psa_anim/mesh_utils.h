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
///\file boundary_marker.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-11-14
///
///\brief

#ifndef PSA_ANIM_TOOLS_PSA_ANIMPY_MESH_UTILS_H
#define PSA_ANIM_TOOLS_PSA_ANIMPY_MESH_UTILS_H

#include <hermes/geometry/bbox.h>
#include <hermes/geometry/plane.h>
#include <hermes/geometry/ray.h>
#include <set>

namespace psa_anim {

bool polyFaceIntersect(const std::vector<hermes::point3> &vertices,
                       const hermes::ray3 &ray);

///
struct FaceMesh {
  std::vector<hermes::point3> vertices;
  std::vector<std::vector<size_t>> faces;

  [[nodiscard]] bool edgeIsSorted(size_t face_id, size_t a, size_t b) const {
    HERMES_ASSERT(face_id < faces.size());
    const auto &face_vertices = faces[face_id];
    int a_index = -1, b_index = -1;
    for (int i = 0; i < face_vertices.size(); ++i) {
      if (face_vertices[i] == a)
        a_index = i;
      if (face_vertices[i] == b)
        b_index = i;
    }
    return b_index == (a_index + 1) % face_vertices.size();
  }

  hermes::point3 faceCenter(size_t face_id) const {
    HERMES_ASSERT(face_id < faces.size());
    const auto &face_vertices = faces[face_id];
    hermes::point3 center;
    for (auto v : face_vertices) {
      center += hermes::vec3(vertices[v]);
    }
    return center / face_vertices.size();
  }
};
/// \note The first id in the pair is always the largest
/// \param faces
/// \param cells
/// \return
std::unordered_map<size_t, std::pair<i64, i64>>
computeFaceTopology(const std::vector<std::vector<size_t>> &faces,
                    const std::vector<std::vector<size_t>> &cells,
                    const std::set<size_t> &region = {});
///
/// \param face_topology
/// \param cells
/// \param cell_id
/// \param ring_size
/// \return
std::vector<size_t> getNeighbours(
    const std::unordered_map<size_t, std::pair<i64, i64>> &face_topology,
    const std::vector<std::vector<size_t>> &cells, size_t cell_id,
    size_t ring_size = 1);

std::vector<size_t>
intersect(const std::unordered_map<size_t, std::pair<i64, i64>> &face_topology,
          const std::vector<std::vector<size_t>> &cells,
          const std::vector<hermes::point3> &centers, size_t cell_id,
          const hermes::bbox3 &region);

///
/// \param face_topology
/// \param cells
/// \param cell_id
/// \param n
/// \return
std::vector<size_t>
getKnn(const std::unordered_map<size_t, std::pair<i64, i64>> &face_topology,
       const std::vector<std::vector<size_t>> &cells,
       const std::vector<hermes::point3> &centers, size_t cell_id, size_t n);

///
/// \param a
/// \param b
/// \return
real_t directionAlignment(const hermes::vec3 &a, const hermes::vec3 &b);
///
/// \param a
/// \param b
/// \return
bool directionIsTheSame(const hermes::vec3 &a, const hermes::vec3 &b);
///
/// \param vertices
/// \param face_vertices
/// \return
hermes::point3 faceCenter(const std::vector<hermes::point3> &vertices,
                          const std::vector<size_t> &face_vertices);
/// \brief Computes cell center (geometric centroid)
/// \param vertices
/// \param faces
/// \param cell_faces
/// \return
hermes::point3 cellCenter(const std::vector<hermes::point3> &vertices,
                          const std::vector<std::vector<size_t>> &faces,
                          const std::vector<size_t> &cell_faces);
///
/// \param vertices
/// \param faces
/// \param cells
/// \return
std::vector<hermes::point3>
cellCenters(const std::vector<hermes::point3> &vertices,
            const std::vector<std::vector<size_t>> &faces,
            const std::vector<std::vector<size_t>> &cells);
/// \brief Computes face's normal using the right-hand rule
/// \param vertices
/// \param face_vertices
/// \return
hermes::vec3 faceNormal(const std::vector<hermes::point3> &vertices,
                        const std::vector<size_t> &face_vertices);
/// \brief Fix faces so boundary faces point outwards their respective cells
/// \param vertices
/// \param cells
/// \param faces
void fixBoundaryFaces(const std::vector<hermes::point3> &vertices,
                      const std::vector<std::vector<size_t>> &cells,
                      std::vector<std::vector<size_t>> &faces);
/// \brief Separate faces of a mesh based on their normals.
/// \note After per face boundary labeling, faces are sorted by boundary index
/// and faces indices are recomputed. \param vertices \param cells \param faces
/// \param boundary_directions
std::vector<std::pair<size_t, size_t>>
markBoundary(const std::vector<hermes::point3> &vertices,
             std::vector<std::vector<size_t>> &cells,
             std::vector<std::vector<size_t>> &faces,
             const std::vector<std::vector<hermes::vec3>> &boundary_directions);
///
/// \param mesh
/// \param directions
/// \return
std::map<std::string, FaceMesh> splitBoundary(
    const FaceMesh &mesh,
    const std::map<std::string, std::vector<hermes::vec3>> &boundary_directions,
    const std::map<std::string, std::vector<std::pair<size_t, size_t>>>
        &boundary_ranges = {});
/// \brief Computes an edge key from two vertex indices
/// \param va
/// \param vb
/// \return
std::pair<size_t, size_t> edgeKey(size_t va, size_t vb);
///
/// \param faces
/// \return
std::map<std::pair<size_t, size_t>, std::pair<i64, i64>>
computeEdgeTopology(const std::vector<std::vector<size_t>> &faces);
///
/// \param v
void flipOrder(std::vector<size_t> &v);
///
FaceMesh genXYAlignedMesh(const FaceMesh &surface,
                          const hermes::point2 &slope_points,
                          float cell_size = 5);
/// \brief Extrudes the surface in the +z direction
/// \note This function do some extra work:
/// \note - Faces are re-ordered to make normals point outwards
/// \param surface
/// \param height
/// \param slope_direction
/// \param slope_point
/// \param boundary_faces_are_triangles
/// \param extrusion_surface_is_inclined
/// \return
FaceMesh extrudeZ(const FaceMesh &surface, real_t height = 100,
                  hermes::vec2 slope_direction = {1, 0},
                  hermes::point2 slope_point = {0, 0},
                  bool extrusion_surface_is_inclined = true,
                  bool generate_boundary_faces = true,
                  bool boundary_faces_are_triangles = true);
/// \brief Loads an obj surface
/// \param path
/// \return
hermes::Result<FaceMesh> loadOBJ(const hermes::Path &path);
///
/// \param mesh
/// \param path
/// \return
bool saveOBJ(const FaceMesh &mesh, const hermes::Path &path,
             bool wireframe = false);
///
/// \param mesh
/// \param path
/// \return
bool saveSTL(const std::map<std::string, FaceMesh> &mesh,
             const hermes::Path &path);
///
/// \param buffer
/// \return
hermes::Result<FaceMesh> parseSTL(const std::string &buffer);
///
/// \param path
/// \return
hermes::Result<FaceMesh> loadSTL(const hermes::Path &path);
/// \brief Loads a stl file containing a surface mesh
/// \param solid
/// \return
hermes::Result<FaceMesh> loadSTL(const hermes::Path &path,
                                 const std::string &solid);

FaceMesh convertToTri(const FaceMesh &input);
/// \brief Merge every two triangles into a quad.
/// \note This function assumes:
/// \note - the triangular mesh was created from a grid where each cell was
/// split into two triangles \note - Each pair of triangles from a respective
/// cell shares the diagonal of that cell \note - All triangles are ordered and
/// their normals point on the same direction \param input \return
FaceMesh convertToQuad(const FaceMesh &input);
/// Projects a grid mesh into another mesh (terrain)
/// \param input
/// \param axis
/// \param origin
/// \param width
/// \param cell_size
/// \return
FaceMesh extractGrid(const FaceMesh &input, const hermes::vec2 &axis,
                     const hermes::point2 &origin, const hermes::vec2 &size,
                     real_t cell_size, bool twoDim = false);
///
/// \param surface
void info(const FaceMesh &surface);

} // namespace psa_anim

#endif // PSA_ANIM_TOOLS_PSA_ANIMPY_MESH_UTILS_H
