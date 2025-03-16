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
///\file terrain_extruder.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-11-30
///
///\brief

#include "mesh_utils.h"
#include <terrain_extruder.h>

namespace psa_anim {

const FaceMesh &
TerrainExtruder::operator[](const std::string &patch_name) const {
  static FaceMesh dummy{};
  if (mesh_.count(patch_name))
    return mesh_.find(patch_name)->second;
  return dummy;
}

void TerrainExtruder::setSlopeDirection(const hermes::point2 &a,
                                        const hermes::point2 &b) {
  slope_direction_ = hermes::normalize(b - a);
  slope_point_ = a;
}

bool TerrainExtruder::load(const hermes::Path &stl_file,
                           const std::string &name, bool align_faces) {
  hermes::Result<FaceMesh> r;
  if (name.empty())
    r = loadSTL(stl_file);
  else
    r = loadSTL(stl_file, name);
  if (!r) {
    HERMES_LOG_ERROR("failed to load stl file: {}", stl_file);
    return false;
  }
  if (align_faces)
    stl_surface_ = genXYAlignedMesh(*r, slope_point_);
  else
    stl_surface_ = *r;
  std::map<std::string, FaceMesh> M;
  M["terrain"] = stl_surface_;
  extrude();
  return true;
}

void TerrainExtruder::save(const hermes::Path &output_path) {
  saveSTL(mesh_, output_path);
}

void TerrainExtruder::extrude(bool top_is_inclined) {
  auto m = extrudeZ(stl_surface_, height, slope_direction_, slope_point_,
                    top_is_inclined, true, true);
  // generate boundary
  auto x = slope_direction_;
  auto y = x.right();
  mesh_ = splitBoundary(m,
                        {{"top", {{0, 0, 1}}},
                         {"outlet", {{x.x, x.y, 0}}},
                         {"inlet", {{-x.x, -x.y, 0}}},
                         {"walls", {{y.x, y.y, 0}, {-y.x, -y.y, 0}}}},
                        {{"terrain", {{0, stl_surface_.faces.size()}}}});
}

} // namespace psa_anim
