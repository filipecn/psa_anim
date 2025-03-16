// Copyright (c) 2022, FilipeCN.
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
///\date 2022-03-11
///
///\brief

#include <hermes/common/arg_parser.h>
#include <hermes/geometry/transform.h>
#include <terrain_extruder.h>

int main(int argc, const char **argv) {
  hermes::Log::addOptions(hermes::logging_options::abbreviate);
  hermes::ArgParser parser("stl2foam",
                           "converts a stl surface into a three-dimensional "
                           "volume for OpenFOAM input.");
  parser.addArgument("-i", "path to the .stl file containing the surface.",
                     true);
  parser.addArgument("-o", "output stl file containing the domain patches.",
                     true);
  parser.addArgument("-p", "input surface name (in stl file).");
  parser.addArgument("--slope-axis-a",
                     "slope main axis begin point in xy plane.");
  parser.addArgument("--slope-axis-b",
                     "slope main axis end point in xy plane.");
  parser.addArgument("--inclined-top", "top follows terrain inclination");
  parser.addArgument("--height", "extrusion height");
  parser.addArgument("--align-faces", "generate an x-aligned mesh");
  HERMES_ASSERT(parser.parse(argc, argv, true));

  // parameters
  hermes::Path path(parser.get<std::string>("-i"));
  hermes::Path output_path(parser.get<std::string>("-o"));
  auto terrain_patch_name = parser.get<std::string>("-p", "terrain");
  auto x_axis_a = parser.getList<real_t>("--slope-axis-a", {0.f, 0.f});
  auto x_axis_b = parser.getList<real_t>("--slope-axis-b", {1.f, 0.f});
  auto top_is_inclined = parser.check("--inclined-top");
  auto height = parser.get<real_t>("--height", 10);
  bool align_faces = parser.check("--align-faces");

  psa_anim::TerrainExtruder extruder;
  extruder.height = height;
  extruder.setSlopeDirection({x_axis_a[0], x_axis_a[1]},
                             {x_axis_b[0], x_axis_b[1]});

  HERMES_ASSERT(extruder.load(path, terrain_patch_name, align_faces));
  extruder.extrude(top_is_inclined);
  extruder.save(output_path);
  HERMES_LOG("...finished stl2foam");
  return 0;
}
