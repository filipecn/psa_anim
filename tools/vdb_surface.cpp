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
///\file vdb_surface.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-03-11
///
///\brief This program generates an iso-surface from each OpenVDB volume in a
/// given folder

#include <hermes/common/arg_parser.h>
#include <hermes/common/file_system.h>
#include <hermes/geometry/transform.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/FastSweeping.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/Statistics.h>
#include <openvdb/tools/VolumeToMesh.h>

int main(int argc, const char **argv) {
  hermes::ArgParser parser("vdb_surface",
                           "extracts a surface mesh from a set of vdb volumes");
  parser.addArgument("-i", "vdb files directory", true);
  parser.addArgument("-o", "output directory", true);
  parser.addArgument("--grid-name", "vdb grid name", true);
  parser.addArgument("--isovalue");
  parser.addArgument("--adaptive");
  parser.addArgument("--smooth");
  parser.parse(argc, argv, true);

  // parameters
  hermes::Path path(parser.get<std::string>("-i"));
  hermes::Path output_path(parser.get<std::string>("-o"));
  auto grid_name = parser.get<std::string>("--grid-name");
  bool adaptive = parser.check("--adaptive");
  bool smooth = parser.check("--smooth");
  auto isovalue = parser.get<real_t>("--isovalue", 0.001f);

  // check output directory
  if (!output_path.isDirectory()) {
    HERMES_LOG_ERROR("invalid output path");
    return -1;
  }

  openvdb::initialize();

  // list all vdb files in the input path
  std::vector<hermes::Path> ls;
  if (path.isFile())
    ls.emplace_back(path);
  else
    ls = hermes::FileSystem::ls(path, hermes::ls_options::files |
                                          hermes::ls_options::sort);

  for (const auto &vdb_file : ls) {
    if (vdb_file.extension() != "vdb")
      continue;

    // read volume grid
    openvdb::io::File file(vdb_file.fullName());
    file.open();
    openvdb::GridBase::Ptr base_grid = file.readGrid(grid_name);
    openvdb::FloatGrid::Ptr grid =
        openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);
    grid->setGridClass(openvdb::GRID_LEVEL_SET);
    file.close();
    if (!grid->activeVoxelCount())
      continue;

    // create level set from avalanche body
    //    auto voxel_size = 0.2;
    //    auto half_width = 1;
    //    openvdb::FloatGrid::Ptr level_set =
    //    openvdb::createLevelSet<openvdb::FloatGrid>(voxel_size, half_width);
    //    for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn();
    //    iter; ++iter)
    //      level_set->fill(iter.getBoundingBox(), 1);

    // create mask grid
    openvdb::FloatGrid::Ptr mask_grid = openvdb::FloatGrid::create();
    mask_grid->denseFill(grid->evalActiveVoxelBoundingBox(), 1);
    for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter;
         ++iter)
      mask_grid->denseFill(iter.getBoundingBox(), 2);
    mask_grid->setGridClass(openvdb::GRID_LEVEL_SET);

    // log
    HERMES_LOG_VARIABLE(vdb_file.fullName());
    // grid->print();

    auto stats = openvdb::v11_0::tools::statistics(grid->cbeginValueOn());
    stats.print();

    auto iv =
        hermes::Numbers::clamp((double)isovalue, stats.mean(), stats.mean());
    HERMES_LOG_VARIABLES(isovalue, stats.max(), stats.mean());
    auto sdf = openvdb::v11_0::tools::fogToSdf(*mask_grid, 1.5);

    openvdb::v11_0::tools::LevelSetFilter<openvdb::FloatGrid> filter(*sdf);
    if (smooth) {
      filter.dilate(5);   // expand the half_width so we can do operations
      filter.gaussian(1); // blur radius in voxels; takes an optional mask as
                          // second parameter
      filter.erode(5);    // reduce the half_width
    }

    // generate mesh
    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec4I> quads;
    std::vector<openvdb::Vec3I> triangles;
    if (adaptive)
      openvdb::v11_0::tools::volumeToMesh(*sdf, points, triangles, quads,
                                          isovalue);
    else
      openvdb::v11_0::tools::volumeToMesh(*sdf, points, quads, isovalue);
    if (points.empty()) {
      HERMES_LOG_WARNING("empty surface");
      continue;
    }

    // store surface
    hermes::Str surf_file_content;
    surf_file_content.appendLine("#v_f3_f4: ", points.size(), " ",
                                 triangles.size(), " ", quads.size());
    for (const auto &point : points)
      surf_file_content.appendLine("v ", point.x(), " ", point.y(), " ",
                                   point.z());
    for (const auto &triangle : triangles)
      surf_file_content.appendLine("f ", triangle.x() + 1, " ",
                                   triangle.y() + 1, " ", triangle.z() + 1);
    for (const auto &quad : quads)
      surf_file_content.appendLine("f ", quad.x() + 1, " ", quad.y() + 1, " ",
                                   quad.z() + 1, " ", quad.w() + 1);

    // construct obj
    auto filename =
        hermes::FileSystem::basename(vdb_file.fullName(), ".vdb") + ".obj";
    HERMES_LOG("output file: [{}]", filename);
    hermes::FileSystem::writeFile(output_path / filename,
                                  surf_file_content.str());
  }
  return 0;
}
