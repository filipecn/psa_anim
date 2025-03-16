/// Copyright (c) 2025, FilipeCN.
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
///\file vdb_edit.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2025-01-10
///
///\brief

#include <hermes/common/arg_parser.h>
#include <hermes/common/debug.h>
#include <hermes/common/file_system.h>
#include <hermes/geometry/plane.h>
#include <hermes/logging/logging.h>
#include <hermes/numeric/numeric.h>
#include <lodepng/lodepng.h>
#include <mesh_utils.h>
#include <openvdb/openvdb.h>
#include <utils.h>

struct Parameters {
  hermes::Path output_path;
  std::vector<hermes::Path> vdb_files;
  bool clip{false};
  hermes::Plane clip_plane;
  hermes::Result<psa_anim::FaceMesh> color_footprint_mesh;
  psa_anim::Color color;
  bool only_missing_frames{false};
  std::vector<unsigned char> image;
  hermes::point2 image_origin{-35.14, -161.9};
  hermes::vec2 image_size{1153, 353};
  hermes::size2 image_pixel_size;

  int read(int argc, const char **argv) {
    hermes::ArgParser parser("vdb_bash", "OpenVDB stats");
    parser.addArgument("-i", "OpenVDB file or directory", true);
    parser.addArgument("-o", "OpenVDB output path");
    parser.addArgument("--mask", "Color mesh");
    parser.addArgument("--clip-plane", "Clip plane");
    parser.addArgument("--image", "color image");
    parser.addArgument("--image-origin", "(x,y)");
    parser.addArgument("--pixel-size", "world size");
    parser.addArgument("--color", "r,g,b");
    parser.addArgument("--only-missing-frames", "process only missing files");
    parser.addArgument("--every-n", "edit only every other nth");
    parser.parse(argc, argv, true);

    output_path = parser.get<std::string>("-o");
    if (!output_path.isDirectory()) {
      HERMES_LOG_ERROR("invalid output path");
      return 0;
    }

    int every_n = parser.get<int>("--every-n", 1);
    only_missing_frames = parser.check("--only-missing-frames");

    hermes::Path path(parser.get<std::string>("-i"));
    if (path.isDirectory()) {
      int k = 0;
      for (const auto &file : hermes::FileSystem::ls(
               path, hermes::ls_options::files | hermes::ls_options::sort))
        if (file.extension() == "vdb") {
          k++;
          if (k % every_n != 0)
            continue;
          if (only_missing_frames &&
              hermes::FileSystem::isFile(output_path / file.name()))
            continue;
          vdb_files.emplace_back(file);
        }
    } else
      vdb_files.emplace_back(path);

    std::string image_str = parser.get<std::string>("--image", "");
    if (!image_str.empty()) {
      float pixel_size = parser.get<float>("--pixel-size", 1.f);
      unsigned error =
          lodepng::decode(image, image_pixel_size.width,
                          image_pixel_size.height, image_str.c_str());
      HERMES_LOG_WARNING("loaded png: {}x{}", image_pixel_size.width,
                         image_pixel_size.height);
      image_size.x = image_pixel_size.width * pixel_size;
      image_size.y = image_pixel_size.height * pixel_size;
      auto origin = parser.getList<float>("--image-origin");
      if (origin.size() == 2) {
        image_origin.x = origin[0];
        image_origin.y = origin[1];
      }
      HERMES_LOG_VARIABLE(image_origin);

      if (error) {
        std::cout << "decoder error " << error << ": "
                  << lodepng_error_text(error) << std::endl;
        return 0;
      }
    }

    auto rgb = parser.getList<int>("--color");
    if (rgb.size() == 3) {
      color.r = rgb[0] / 255.;
      color.g = rgb[1] / 255.;
      color.b = rgb[2] / 255.;
    }

    auto clip_values = parser.getList<float>("--clip-plane");
    if (clip_values.size() == 6) {
      clip_plane =
          hermes::Plane({clip_values[0], clip_values[1], clip_values[2]},
                        {clip_values[3], clip_values[4], clip_values[5]});
      clip = true;
    }

    std::string color_footprint_str = parser.get<std::string>("--mask", "");
    if (!color_footprint_str.empty()) {
      color_footprint_mesh = psa_anim::loadOBJ(color_footprint_str);
    }

    return 1;
  }
} args;

psa_anim::Color imageColor(const openvdb::Vec3f &p) {
  if (args.image.empty())
    return {1.0, 1.0, 0.0, 1.0};

  hermes::vec2 uv = hermes::point2(p.x(), p.y()) - args.image_origin;
  uv /= args.image_size;
  uv.x = hermes::Numbers::clamp(uv.x, 0.f, 1.f);
  uv.y = hermes::Numbers::clamp(uv.y, 0.f, 1.f);
  uv.y = 1 - uv.y;

  hermes::index2 pixel(uv.x * args.image_pixel_size.width,
                       uv.y * args.image_pixel_size.height);
  pixel.i =
      hermes::Numbers::clamp(pixel.i, 0, (int)args.image_pixel_size.width - 1);
  pixel.j =
      hermes::Numbers::clamp(pixel.j, 0, (int)args.image_pixel_size.height - 1);

  return {
      args.image[pixel.j * args.image_pixel_size.width * 4 + pixel.i * 4 + 0] /
          255.f,
      args.image[pixel.j * args.image_pixel_size.width * 4 + pixel.i * 4 + 1] /
          255.f,
      args.image[pixel.j * args.image_pixel_size.width * 4 + pixel.i * 4 + 2] /
          255.f,
      args.image[pixel.j * args.image_pixel_size.width * 4 + pixel.i * 4 + 3] /
          255.f};
}

bool footprintColor(const openvdb::Vec3f &p) {
  if (args.color_footprint_mesh->faces.empty())
    return false;
  for (const auto &face : args.color_footprint_mesh->faces) {
    size_t face_size = face.size();
    bool side = false;
    bool is_inside = true;
    for (size_t i = 0; i < face_size && is_inside; ++i) {
      auto a =
          args.color_footprint_mesh->vertices[face[(i + 0) % face_size]].xy();
      auto b =
          args.color_footprint_mesh->vertices[face[(i + 1) % face_size]].xy();
      bool this_side =
          hermes::cross(a - b, hermes::point2(p.x(), p.y()) - b) >= 0;
      if (i == 0)
        side = this_side;
      else if (side != this_side)
        is_inside = false;
    }
    if (is_inside)
      return true;
  }
  return false;
}

void process(const hermes::Path &path) {
  HERMES_LOG("process: {}", path);
  // load vdb
  openvdb::io::File file(path.fullName());
  file.open();

  std::set<std::string> grid_names;
  for (openvdb::io::File::NameIterator name_iter = file.beginName();
       name_iter != file.endName(); ++name_iter)
    grid_names.insert(name_iter.gridName());

  HERMES_ASSERT(grid_names.count("density"));

  openvdb::GridBase::Ptr base_grid = file.readGrid("density");
  openvdb::GridBase::Ptr color_base_grid;

  openvdb::FloatGrid::Ptr grid =
      openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);

  openvdb::Vec3fGrid::Ptr color_vdb;
  if (grid_names.count("color")) {
    color_base_grid = file.readGrid("color");
    color_vdb = openvdb::gridPtrCast<openvdb::Vec3fGrid>(color_base_grid);
  } else {
    color_vdb = openvdb::Vec3fGrid::create();
    color_vdb->setTransform(grid->transformPtr());
    color_vdb->setName("color");
  }

  file.close();

  openvdb::Vec3fGrid::Accessor a = color_vdb->getAccessor();

  for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter.test();
       ++iter) {
    openvdb::Vec3f wp = grid->transform().indexToWorld(iter.getCoord());
    if (footprintColor(wp))
      a.setValue(iter.getCoord(), {args.color.r, args.color.g, args.color.b});
    if (!args.image.empty()) {
      auto color = imageColor(wp);
      a.setValue(iter.getCoord(), {color.r, color.g, color.b});
    }
    if (args.clip &&
        !args.clip_plane.onNormalSide(hermes::point3(wp.x(), wp.y(), wp.z()))) {
      iter.setValue(0);
      iter.setValueOff();
    }
  }

  auto output_file = args.output_path / path.name();
  HERMES_LOG("writting {}", output_file.fullName());

  openvdb::io::File(output_file.fullName()).write({grid, color_vdb});
}

int main(int argc, const char **argv) {
  openvdb::initialize();

  if (!args.read(argc, argv))
    return -1;

  // process
  for (const auto &path : args.vdb_files) {
    process(path);
  }

  return 0;
}
