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
///\file block_mesh.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-11-16
///
///\brief

#ifndef PSA_ANIM_TOOLS_PSA_ANIMPY_BLOCK_MESH_H
#define PSA_ANIM_TOOLS_PSA_ANIMPY_BLOCK_MESH_H

#include <hermes/common/file_system.h>
#include <hermes/geometry/point.h>
#include <mesh_utils.h>

namespace psa_anim {

// *********************************************************************************************************************
//                                                                                               BlockMeshDescription
// *********************************************************************************************************************
/// \brief OpenFOAM's block mesh description
class BlockMeshDescription {
public:
  struct Face {
    std::vector<size_t> vertices;
  };
  struct Edge {
    std::string type;
    size_t a;
    size_t b;
    std::vector<hermes::point3> interpolation_points;
  };
  struct Block {
    std::vector<size_t> hex;
    hermes::vec3 grading{1, 1, 1};
    hermes::size3 resolution{1, 1, 1};

    ///
    /// \param i
    /// \return
    std::vector<size_t> getFace(size_t i) const {
      // blocks follow the vertex sequence
      //                       f3
      //                 7 _______  6
      //                /  f5    /
      //            4  _______5/  |  f1
      //              |       |   |
      // x3   x2  f0  |  f2   |   / 2
      //  | /         |       | /
      //   -- x1     0  ------ 1
      //                 f4
      static std::vector<std::vector<size_t>> indices = {
          {0, 4, 7, 3}, // f0
          {1, 2, 6, 5}, // f1
          {0, 1, 5, 4}, // f2
          {2, 3, 7, 6}, // f3
          {0, 3, 2, 1}, // f4
          {4, 5, 6, 7}, // f5
      };
      HERMES_ASSERT(i < 6);
      HERMES_ASSERT(hex.size() == 8);
      std::vector<size_t> face_vertices(4);
      for (size_t j = 0; j < 4; ++j)
        face_vertices[j] = hex[indices[i][j]];
      return face_vertices;
    }
  };
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  ///
  /// \param v
  void pushVertex(const hermes::point3 &v) { vertices_.emplace_back(v); }
  ///
  /// \param edge
  void pushEdge(const Edge &edge) { edges_.emplace_back(edge); }
  ///
  /// \param block
  void pushBlock(const Block &block) { blocks_.emplace_back(block); }
  ///
  /// \param name
  /// \param type
  void setFaceType(const std::string &name, const std::string &type) {
    boundary_types_[name] = type;
  }
  /// \param name
  /// \param face
  void addFaceToBoundary(const std::string &name,
                         const std::vector<size_t> &face) {
    boundary_faces_[name].emplace_back(face);
  }
  ///
  /// \param path
  void save(const hermes::Path &path) const {
    hermes::Str content;
    content.appendLine(
        "/*--------------------------------*- C++ "
        "-*----------------------------------*\\\n"
        "  =========                 |\n"
        "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD "
        "Toolbox\n"
        "   \\\\    /   O peration     | Website:  https://openfoam.org\n"
        "    \\\\  /    A nd           | Version:  6\n"
        "     \\\\/     M anipulation  |\n"
        "\\*-------------------------------------------------------------------"
        "--------*/\n"
        "FoamFile\n"
        "{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        "    class       dictionary;\n"
        "    object      blockMeshDict;\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
        "* * * * //\n\n");

    content.appendLine("scale ", convert_to_meters_scale, ";");

    // write vertices
    {
      content.appendLine("vertices");
      content.appendLine("(");
      for (size_t i = 0; i < vertices_.size(); ++i)
        content.appendLine("\t(", vertices_[i].x, " ", vertices_[i].y, " ",
                           vertices_[i].z, ")\t// vertex number ", i);
      content.appendLine(");\n\n");
    }

    // write blocks
    {
      content.appendLine("blocks");
      content.appendLine("(");
      for (size_t i = 0; i < blocks_.size(); ++i) {
        const auto &block = blocks_[i];
        content.append("\thex (", hermes::Str::join(block.hex, " "), ") ");
        content.append("(", block.resolution.width, " ",
                       block.resolution.height, " ", block.resolution.depth,
                       ") ");
        content.append("simpleGrading (", block.grading.x, " ", block.grading.y,
                       " ", block.grading.z, ")");
        content.appendLine("\t//BLOCK ", i);
      }
      content.appendLine(");\n\n");
    }

    // write boundary
    {
      content.appendLine("boundary");
      content.appendLine("(");
      for (const auto &patch : boundary_faces_) {
        content.appendLine("\t", patch.first);
        content.appendLine("\t{");
        std::string type = "patch";
        auto patch_it = boundary_types_.find(patch.first);
        if (patch_it != boundary_types_.end())
          type = patch_it->second;
        content.appendLine("\t\ttype ", type, ";");
        content.appendLine("\t\tfaces");
        content.appendLine("\t\t(");
        for (const auto &face : patch.second)
          content.appendLine("\t\t(", hermes::Str::join(face, " "), ")");
        content.appendLine("\t\t);");
        content.appendLine("\t}");
      }
      content.appendLine(");\n");
    }

    // write file
    hermes::FileSystem::writeFile(path, content.str());
  }

  void saveOBJ(const hermes::Path &path) const {
    hermes::Str s;
    for (const auto &v : vertices_)
      s.appendLine("v ", v.x, " ", v.y, " ", v.z);
    for (const auto &block : blocks_)
      for (size_t i = 0; i < 6; ++i) {
        auto face = block.getFace(i);
        for (auto &f : face)
          f += 1;
        s.appendLine("f ", hermes::Str::join(face, " "));
      }
    hermes::FileSystem::writeFile(path, s.str());
  }

  void info() const {
    HERMES_LOG("Block Mesh Descriptor");
    HERMES_LOG_VARIABLE(offset);
    HERMES_LOG_VARIABLE(vertices_.size());
    HERMES_LOG_VARIABLE(edges_.size());
    HERMES_LOG_VARIABLE(blocks_.size());
    for (const auto &boundary : boundary_faces_)
      HERMES_LOG("#boundary[{}] = {}", boundary.first, boundary.second.size());
  }

  const hermes::point3 &vertexAt(size_t i) const {
    HERMES_ASSERT(i < vertices_.size());
    return vertices_[i];
  }

  // *******************************************************************************************************************
  //                                                                                                    PUBLIC FIELDS
  // *******************************************************************************************************************
  real_t convert_to_meters_scale{1.0};
  hermes::vec3 offset;

private:
  std::vector<hermes::point3> vertices_;
  std::vector<Edge> edges_;
  std::vector<Block> blocks_;
  std::unordered_map<std::string, std::vector<std::vector<size_t>>>
      boundary_faces_;
  std::unordered_map<std::string, std::string> boundary_types_;
};

struct BlockMeshConfig {
  bool face_matching = true;
  hermes::vec3 grading{1, 1, 1};
  hermes::size3 resolution{1, 1, 1};
};

/// \brief Generates a block mesh from a terrain surface
/// \note This function extrudes the faces of a surface and produce a block for
/// each
/// \note Input boundaries MUST contain the patches "top" and "bottom" (or
/// "terrain")
/// \param vertices
/// \param faces quad mesh
/// \param boundaries direction of boundaries
/// \param top_is_inclined
/// \return
BlockMeshDescription
blockMesh(const std::vector<hermes::point3> &vertices,
          const std::vector<std::vector<size_t>> &faces,
          std::map<std::string, std::vector<hermes::vec3>> boundaries,
          real_t height = 100, BlockMeshConfig config = {},
          hermes::vec2 slope_direction = {1, 0},
          hermes::point2 slope_point = {0, 0}, bool top_is_inclined = true);

} // namespace psa_anim

#endif // PSA_ANIM_TOOLS_PSA_ANIMPY_BLOCK_MESH_H
