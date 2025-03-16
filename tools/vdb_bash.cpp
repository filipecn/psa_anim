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
///\file foam2vdb.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-03-11
///
///\brief

#include <hermes/common/arg_parser.h>
#include <hermes/common/debug.h>
#include <hermes/common/file_system.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Statistics.h>

std::pair<openvdb::v11_0::math::Stats, openvdb::v11_0::math::Histogram>
process(const hermes::Path &path) {
  openvdb::io::File file(path.fullName());
  file.open();

  std::vector<std::string> grid_names;
  for (openvdb::io::File::NameIterator name_iter = file.beginName();
       name_iter != file.endName(); ++name_iter)
    grid_names.emplace_back(name_iter.gridName());

  auto grid_name = grid_names[0];

  openvdb::GridBase::Ptr base_grid = file.readGrid(grid_name);

  file.close();

  openvdb::FloatGrid::Ptr grid =
      openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);

  auto stats = openvdb::v11_0::tools::statistics(grid->cbeginValueOn());

  return {stats, openvdb::v11_0::tools::histogram(
                     grid->cbeginValueAll(), stats.min(), stats.max(), 20)};
}

int main(int argc, const char **argv) {
  hermes::ArgParser parser("vdb_bash", "OpenVDB stats");
  parser.addArgument("-i", "OpenVDB file or directory", true);
  parser.parse(argc, argv, true);
  // parameters
  hermes::Path path(parser.get<std::string>("-i"));
  std::vector<hermes::Path> vdb_files;
  if (path.isDirectory()) {
    for (const auto &file : hermes::FileSystem::ls(
             path, hermes::ls_options::files | hermes::ls_options::sort))
      if (file.extension() == "vdb")
        vdb_files.emplace_back(file);
  } else
    vdb_files.emplace_back(path);

  openvdb::initialize();

  size_t highest_mean_index = 0;
  size_t highest_density_index = 0;
  double highest_density = 0;
  double highest_mean = 0;

  for (size_t i = 0; i < vdb_files.size(); ++i) {
    HERMES_LOG("file[{}]: {}", i, vdb_files[i].fullName());
    auto p = process(vdb_files[i]);

    if (highest_density < p.first.max()) {
      highest_density = p.first.max();
      highest_density_index = i;
    }

    if (highest_mean < p.first.mean()) {
      highest_mean = p.first.mean();
      highest_mean_index = i;
    }
  }

  auto p = process(vdb_files[highest_mean_index]);
  p.first.print();
  p.second.print();

  p = process(vdb_files[highest_density_index]);
  p.first.print();
  p.second.print();

  return 0;
}
