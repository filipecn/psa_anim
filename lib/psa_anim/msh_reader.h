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
///\file msh_reader.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-11-11
///
///\brief

#ifndef PSA_ANIM_TOOLS_PSA_ANIMPY_MSH_READER_H
#define PSA_ANIM_TOOLS_PSA_ANIMPY_MSH_READER_H

#include <hermes/common/debug.h>
#include <hermes/common/file_system.h>
#include <hermes/common/parsers.h>
#include <hermes/data_structures/multi_hash_map.h>
#include <hermes/geometry/point.h>
#include <set>

namespace psa_anim {

// *********************************************************************************************************************
//                                                                                                            MshMesh
// *********************************************************************************************************************
class MshMesh {
public:
  struct Node {
    size_t id;
    hermes::point3 position;
  };
  struct Element {
    size_t id;
    size_t type;
    std::vector<size_t> tags;
    std::vector<size_t> nodes;
  };
  // *******************************************************************************************************************
  //                                                                                                   STATIC METHODS
  // *******************************************************************************************************************
  /// \note Based on http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
  /// \param element_type
  /// \return
  static size_t elementSize(size_t element_type) {
    static size_t geometrical_element_size[32] = {
        0,  // [0]
        2,  // [1]:  2-node line.
        3,  // [2]:  3-node triangle.
        4,  // [3]:  4-node quadrangle.
        4,  // [4]:  4-node tetrahedron.
        8,  // [5]:  8-node hexahedron.
        6,  // [6]:  6-node prism.
        5,  // [7]:  5-node pyramid.
        3,  // [8]:  3-node second order line.
        6,  // [9]:  6-node second order triangle.
        9,  // [10]: 9-node second order quadrangle.
        10, // [11]: 10-node second order tetrahedron.
        27, // [12]: 27-node second order hexahedron.
        18, // [13]: 18-node second order prism.
        14, // [14]: 14-node second order pyramid.
        1,  // [15]: 1-node point.
        8,  // [16]: 8-node second order quadrangle.
        20, // [17]: 20-node second order hexahedron.
        15, // [18]: 15-node second order prism.
        13, // [19]: 13-node second order pyramid.
        9,  // [20]: 9-node third order incomplete triangle.
        10, // [21]: 10-node third order triangle.
        12, // [22]: 12-node fourth order incomplete triangle.
        15, // [23]: 15-node fourth order triangle.
        15, // [24]: 15-node fifth order incomplete triangle.
        21, // [25]: 21-node fifth order complete triangle.
        4,  // [26]: 4-node third order edge.
        5,  // [27]: 5-node fourth order edge.
        6,  // [28]: 6-node fifth order edge.
        20, // [29]: 20-node third order tetrahedron.
        35, // [30]: 35-node fourth order tetrahedron.
        56  // [31]: 56-node fifth order tetrahedron.
    };
    if (element_type > 31)
      return 0;
    return geometrical_element_size[element_type];
  }
  /// \note Based on
  /// http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_65.php#SEC65
  /// \param element_type
  /// \return
  static std::vector<std::vector<size_t>> subElements(size_t element_type) {
    std::vector<std::vector<std::vector<size_t>>> ordering = {
        {{}},           // [0]
        {{0, 1}},       // [1]:  2-node line.
        {{0, 1, 2}},    // [2]:  3-node triangle.
        {{0, 1, 2, 3}}, // [3]:  4-node quadrangle.
        {{0, 2, 1},
         {1, 2, 3},
         {0, 3, 2},
         {0, 1, 3}}, // [4]:  4-node tetrahedron.
                     //        {},       // [5]:  8-node hexahedron.
                     //        {},       // [6]:  6-node prism.
                     //        {},       // [7]:  5-node pyramid.
        //        {},       // [8]:  3-node second order line (2 nodes
        //        associated with the vertices and 1 with the edge).
        //        {},       // [9]:  6-node second order triangle (3 nodes
        //        associated with the vertices and 3 with the edges).
        //        {},       // [10]: 9-node second order quadrangle (4 nodes
        //        associated with the vertices, 4 with the edges and 1 with the
        //        face).
        //        {},       // [11]: 10-node second order tetrahedron (4 nodes
        //        associated with the vertices and 6 with the edges).
        //        {},       // [12]: 27-node second order hexahedron (8 nodes
        //        associated with the vertices, 12 with the edges, 6 with the
        //        faces and 1 with the volume).
        //        {},       // [13]: 18-node second order prism (6 nodes
        //        associated with the vertices, 9 with the edges and 3 with the
        //        quadrangular faces).
        //        {},       // [14]: 14-node second order pyramid (5 nodes
        //        associated with the vertices, 8 with the edges and 1 with the
        //        quadrangular face).
        //        {},       // [15]: 1-node point.
        //        {},       // [16]: 8-node second order quadrangle (4 nodes
        //        associated with the vertices and 4 with the edges).
        //        {},       // [17]: 20-node second order hexahedron (8 nodes
        //        associated with the vertices and 12 with the edges).
        //        {},       // [18]: 15-node second order prism (6 nodes
        //        associated with the vertices and 9 with the edges).
        //        {},       // [19]: 13-node second order pyramid (5 nodes
        //        associated with the vertices and 8 with the edges).
        //        {},       // [20]: 9-node third order incomplete triangle (3
        //        nodes associated with the vertices, 6 with the edges)
        //        {},       // [21]: 10-node third order triangle (3 nodes
        //        associated with the vertices, 6 with the edges, 1 with the
        //        face)
        //        {},       // [22]: 12-node fourth order incomplete triangle (3
        //        nodes associated with the vertices, 9 with the edges)
        //        {},       // [23]: 15-node fourth order triangle (3 nodes
        //        associated with the vertices, 9 with the edges, 3 with the
        //        face)
        //        {},       // [24]: 15-node fifth order incomplete triangle (3
        //        nodes associated with the vertices, 12 with the edges)
        //        {},       // [25]: 21-node fifth order complete triangle (3
        //        nodes associated with the vertices, 12 with the edges, 6 with
        //        the face)
        //        {},       // [26]: 4-node third order edge (2 nodes associated
        //        with the vertices, 2 internal to the edge)
        //        {},       // [27]: 5-node fourth order edge (2 nodes
        //        associated with the vertices, 3 internal to the edge)
        //        {},       // [28]: 6-node fifth order edge (2 nodes associated
        //        with the vertices, 4 internal to the edge)
        //        {},       // [29]: 20-node third order tetrahedron (4 nodes
        //        associated with the vertices, 12 with the edges, 4 with the
        //        faces)
        //        {},       // [30]: 35-node fourth order tetrahedron (4 nodes
        //        associated with the vertices, 18 with the edges, 12 with the
        //        faces, 1 in the volume)
        //        {}        // [31]: 56-node fifth order tetrahedron (4 nodes
        //        associated with the vertices, 24 with the edges, 24 with the
        //        faces, 4 in the volume)
    };
    if (element_type > 4)
      HERMES_NOT_IMPLEMENTED;
    if (element_type > 32)
      return {};
    return ordering[element_type];
  }
  // *******************************************************************************************************************
  //                                                                                                 FRIEND FUNCTIONS
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  //                                                                                                       assignment
  // *******************************************************************************************************************
  //                                                                                                        OPERATORS
  // *******************************************************************************************************************
  //                                                                                                       assignment
  //                                                                                                       arithmetic
  //                                                                                                          boolean
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  /// \note the mesh is appended to the input vectors, i.e, input data is not
  /// cleared
  /// \param vertices
  /// \param faces
  /// \param cells
  void retrieveMesh(std::vector<hermes::point3> &vertices,
                    std::vector<std::vector<size_t>> &faces,
                    std::vector<std::vector<size_t>> &cells) const {
    // map of vertices: msh index -> vertices index
    std::unordered_map<size_t, size_t> msh_node_vertex_index_map;
    // map of faces: face vertices -> face index
    hermes::MultiHashMap<size_t, size_t> face_vertices_to_face_index;
    // iterate over elements to construct & populate the mesh
    for (const auto &element : elements_) {
      std::vector<size_t> cell_faces;
      const auto &sub_elements = subElements(element.type);
      // for each face of the element generate a face
      // also map the vertices and face indices
      for (const auto &sub_element : sub_elements) {
        std::vector<size_t> face_nodes;
        for (size_t face_node_index : sub_element) {
          HERMES_ASSERT(face_node_index != element.nodes.size());
          auto msh_vertex_index = element.nodes[face_node_index];
          auto node_index_it = node_id_node_index_map_.find(msh_vertex_index);
          HERMES_ASSERT(node_index_it != node_id_node_index_map_.end());
          msh_vertex_index = node_index_it->second;
          // register the vertex if this is a new one
          if (msh_node_vertex_index_map.find(msh_vertex_index) ==
              msh_node_vertex_index_map.end()) {
            msh_node_vertex_index_map[msh_vertex_index] = vertices.size();
            HERMES_ASSERT(msh_vertex_index < nodes_.size());
            vertices.emplace_back(nodes_[msh_vertex_index].position);
          }
          // the face node is converted to the vertices index
          face_nodes.emplace_back(msh_node_vertex_index_map[msh_vertex_index]);
        }
        // now we have a list of vertices that compose the current face
        // lets register the face
        // The face is uniquely labeled by these vertices sorted
        auto face_nodes_key = face_nodes;
        std::sort(face_nodes_key.begin(), face_nodes_key.end());
        if (!face_vertices_to_face_index.contains(face_nodes_key)) {
          face_vertices_to_face_index.insert(face_nodes_key, faces.size());
          // add the face
          faces.emplace_back(face_nodes);
        }
        auto face_id_r = face_vertices_to_face_index.get(face_nodes_key);
        HERMES_ASSERT(face_id_r);
        // now that we have registered the vertices of the face, and the face
        // itself we need to add the face to the cell
        cell_faces.emplace_back(*face_id_r);
      }
      cells.emplace_back(cell_faces);
    }
  }
  ///
  /// \return
  bool isTetMesh() const {
    return std::all_of(elements_.begin(), elements_.end(),
                       [](const auto &element) -> bool {
                         if (element.type != 4)
                           return false;
                         return true;
                       });
  }
  ///
  /// \param id
  /// \param position
  void addNode(size_t id, const hermes::point3 &position) {
    node_id_node_index_map_[id] = nodes_.size();
    nodes_.push_back({.id = id, .position = position});
  }
  ///
  /// \param element
  void addElement(const Element &element) { elements_.emplace_back(element); }
  ///
  /// \param id
  /// \param type
  /// \param tags
  /// \param element_nodes
  void addElement(size_t id, size_t type, const std::vector<size_t> &tags,
                  const std::vector<size_t> &element_nodes) {
    elements_.push_back(
        {.id = id, .type = type, .tags = tags, .nodes = element_nodes});
  }
  // *******************************************************************************************************************
  //                                                                                                    PUBLIC FIELDS
  // *******************************************************************************************************************
private:
  std::unordered_map<size_t, size_t> node_id_node_index_map_;
  std::vector<Node> nodes_;
  std::vector<Element> elements_;
};

// *********************************************************************************************************************
//                                                                                                            MshMesh
// *********************************************************************************************************************
class MshReader {
public:
  // *******************************************************************************************************************
  //                                                                                                   STATIC METHODS
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  //                                                                                                 FRIEND FUNCTIONS
  // *******************************************************************************************************************
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  MshReader() {
    msh_parser_.setBlankCharacters(" \t\n");
    msh_parser_.pushBlockDelimiters("\\$MeshFormat", "\\$EndMeshFormat");
    msh_parser_.pushBlockDelimiters("\\$Nodes", "\\$EndNodes");
    msh_parser_.pushBlockDelimiters("\\$Elements", "\\$EndElements");
    msh_parser_.pushTokenPattern("integer", hermes::Str::regex::integer_number);
    msh_parser_.pushTokenPattern("real",
                                 hermes::Str::regex::floating_point_number);
  }
  //                                                                                                       assignment
  // *******************************************************************************************************************
  //                                                                                                        OPERATORS
  // *******************************************************************************************************************
  //                                                                                                       assignment
  //                                                                                                       arithmetic
  //                                                                                                          boolean
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  hermes::Result<MshMesh> loadFromFile(const hermes::Path &path) {
    HERMES_LOG("reading {}", path);
    MshMesh mesh;
#define EXPECT_PREFIX(P)                                                       \
  if (!hermes::Str::isPrefix(P, nextLine())) {                                 \
    HERMES_LOG("expected {} at line [{}]: {}", P, current_line_number,         \
               lines[current_line_number]);                                    \
    return hermes::Result<MshMesh>::error(HeResult::INVALID_INPUT);            \
  }

#define EXPECT(C)                                                              \
  if (!(C)) {                                                                  \
    HERMES_LOG_ERROR("expected {} at line [{}]: {}", #C, current_line_number,  \
                     lines[current_line_number]);                              \
    return hermes::Result<MshMesh>::error(HeResult::INVALID_INPUT);            \
  }

#define READ_NUMBER_OF_ENTRIES(N)                                              \
  auto n_t = hermes::Str::strip(nextLine(), " \n\t");                          \
  EXPECT(hermes::Str::isInteger(n_t))                                          \
  auto N = std::stoull(n_t)

    auto lines = hermes::FileSystem::readLines(path);
    auto number_of_lines = lines.size();
    size_t current_line_number = 0;

    auto currentLine = [&]() -> const std::string & {
      static const std::string empty_line;
      if (current_line_number < number_of_lines)
        return lines[current_line_number];
      return empty_line;
    };

    auto nextLine = [&]() -> const std::string & {
      static const std::string empty_line;
      current_line_number++;
      if (current_line_number < number_of_lines)
        return lines[current_line_number];
      return empty_line;
    };

    while (current_line_number < lines.size()) {
      if (hermes::Str::isPrefix("$MeshFormat", currentLine())) {
        // read mesh format
        HERMES_LOG(nextLine().c_str());
        // now we should find the end of the block
        EXPECT_PREFIX("$EndMeshFormat")
        current_line_number++;
      }
      if (hermes::Str::isPrefix("$Nodes", currentLine())) {
        READ_NUMBER_OF_ENTRIES(number_of_nodes);
        for (size_t i = 0; i < number_of_nodes; ++i) {
          const auto &line = nextLine();
          // here we expect a point id x y z
          auto text_pieces = hermes::Str::split(hermes::Str::strip(line));
          EXPECT(text_pieces.size() == 4);
          EXPECT(hermes::Str::isInteger(text_pieces[0]));
          EXPECT(hermes::Str::isNumber(text_pieces[1]));
          EXPECT(hermes::Str::isNumber(text_pieces[2]));
          EXPECT(hermes::Str::isNumber(text_pieces[3]));
          auto node_id = std::stoull(text_pieces[0]);
          EXPECT(node_id <= number_of_nodes);
          mesh.addNode(node_id,
                       {std::stof(text_pieces[1]), std::stof(text_pieces[2]),
                        std::stof(text_pieces[3])});
        }
        EXPECT_PREFIX("$EndNodes")
        HERMES_LOG("... {} nodes loaded ...", number_of_nodes);
        current_line_number++;
      }
      if (hermes::Str::isPrefix("$Elements", currentLine())) {
        READ_NUMBER_OF_ENTRIES(number_of_elements);
        for (size_t i = 0; i < number_of_elements; ++i) {
          const auto &line = nextLine();
          auto text_pieces = hermes::Str::split(hermes::Str::strip(line));
          for (const auto &text_piece : text_pieces)
            EXPECT(hermes::Str::isInteger(text_piece));
          // elm-number elm-type number-of-tags < tag > … node-number-list
          MshMesh::Element element;
          EXPECT(text_pieces.size() > 3);
          element.id = std::stoull(text_pieces[0]);
          element.type = std::stoull(text_pieces[1]);
          size_t element_size = MshMesh::elementSize(element.type);
          size_t n_tags = std::stoull(text_pieces[2]);
          EXPECT(text_pieces.size() == 3 + n_tags + element_size);
          for (size_t j = 0; j < n_tags; ++j)
            element.tags.emplace_back(std::stoull(text_pieces[3 + j]));
          for (size_t j = 0; j < element_size; ++j)
            element.nodes.emplace_back(
                std::stoull(text_pieces[3 + n_tags + j]));
          mesh.addElement(element);
        }
        EXPECT_PREFIX("$EndElements")
        HERMES_LOG("... {} elements loaded ...", number_of_elements);
        current_line_number++;
      }
      current_line_number++;
    }

#undef EXPECT_PREFIX
#undef EXPECT
#undef READ_NUMBER_OF_ENTRIES

    return hermes::Result<MshMesh>(mesh);
  }
  // *******************************************************************************************************************
  //                                                                                                    PUBLIC FIELDS
  // *******************************************************************************************************************
private:
  hermes::StringParser msh_parser_;
};

} // namespace psa_anim

#endif // PSA_ANIM_TOOLS_PSA_ANIMPY_MSH_READER_H
