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
///\file stl2foam.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-10-31
///
///\brief

#ifndef PSA_ANIM_TOOLS_STL2FOAM_STL2FOAM_H
#define PSA_ANIM_TOOLS_STL2FOAM_STL2FOAM_H

#include <hermes/common/file_system.h>
#include <hermes/geometry/bbox.h>
#include <hermes/geometry/point.h>
#include <mesh_utils.h>
#include <queue>

namespace psa_anim {

// *********************************************************************************************************************
//                                                                                                    TerrainExtruder
// *********************************************************************************************************************
/// \brief Tool for extruding a triangle mesh.
/// \note The extrusion is made towards the up direction.
/// \note 5 patches are generated: terrain, top, walls, inlet and outlet.
/// \note The faces facing the down direction are assigned to the patch named
/// terrain.
/// \note The faces facing the up direction are assigned to the patch named top.
/// \note The faces facing the left direction are assigned to the patch named
/// inlet.
/// \note The faces facing the right direction are assigned to the patch named
/// outlet.
/// \note The remaining faces are assigned to the patch named walls.
class TerrainExtruder {
public:
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  TerrainExtruder() = default;
  //                                                                                                       assignment
  // *******************************************************************************************************************
  //                                                                                                        OPERATORS
  // *******************************************************************************************************************
  ///
  /// \param patch_name
  /// \return
  const FaceMesh &operator[](const std::string &patch_name) const;
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  /// Sets the main direction in XY plane
  void setSlopeDirection(const hermes::point2 &a, const hermes::point2 &b);
  /// \brief Loads a stl file containing a surface mesh
  /// \param stl_file directory or file
  /// \return true on success
  bool load(const hermes::Path &stl_file, const std::string &name = "",
            bool align_faces = false);
  /// \brief Stores the mesh
  /// \param output_path
  void save(const hermes::Path &output_path);
  /// \brief Extrudes the stl_surface mesh
  void extrude(bool top_is_inclined = true);
  // *******************************************************************************************************************
  //                                                                                                    PUBLIC FIELDS
  // *******************************************************************************************************************
  real_t height = 100;

private:
  ///
  std::map<std::string, FaceMesh> mesh_;
  ///
  FaceMesh stl_surface_;
  // up direction
  const hermes::vec3 up_direction_{0, 0, 1};
  // medial axis info
  hermes::vec2 slope_direction_{1, 0};
  hermes::point2 slope_point_{};
};

} // namespace psa_anim

#endif // PSA_ANIM_TOOLS_STL2FOAM_STL2FOAM_H
