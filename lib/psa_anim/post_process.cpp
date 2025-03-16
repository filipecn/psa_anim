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
///\file post_process.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-12-07
///
///\brief

#include <post_process.h>

namespace psa_anim {

double w(const hermes::point3 &p, const hermes::point3 &q) {
  return hermes::distance2(p, q);
}

std::vector<double> smooth(const PolyMesh *mesh,
                           const std::vector<double> &field_values,
                           double isovalue, size_t size) {
  std::vector<double> s_field = field_values;
  if (!mesh)
    return std::move(s_field);
  if (mesh->cells.size() != field_values.size()) {
    HERMES_LOG_ERROR("mesh cells count differs from cell values count!");
    return std::move(s_field);
  }
  HERMES_ASSERT(!mesh->cell_centers.empty())

  auto topology = computeFaceTopology(mesh->faces, mesh->cells);
  // here we detect "interior" cells based on the isovalue
  std::vector<bool> is_inside(field_values.size(), false);
  for (size_t i = 0; i < is_inside.size(); ++i)
    if (field_values[i] >= isovalue)
      is_inside[i] = true;
  // now we smooth the field into every boundary cell
  // a boundary cell is any cell with sufficient field value neighbouring an
  // inside cell
  for (size_t cell_id = 0; cell_id < mesh->cells.size(); ++cell_id) {
    //    if (is_inside[cell_id])
    //      continue;
    // retrieve neighbours
    auto neighbours = getNeighbours(topology, mesh->cells, cell_id, size);
    std::vector<size_t> inside_neighbours;
    for (auto neighbour : neighbours) {
      HERMES_ASSERT(neighbour < field_values.size());
      if (is_inside[neighbour])
        inside_neighbours.emplace_back(neighbour);
    }
    if (!inside_neighbours.empty()) {
      // method 0: mean
      {
        double value_sum = field_values[cell_id];
        for (auto neighbour : inside_neighbours) {
          HERMES_ASSERT(neighbour < field_values.size());
          value_sum += field_values[neighbour];
        }
        s_field[cell_id] = value_sum / (inside_neighbours.size() + 1);
      }
      continue;
      // method 1: laplacian diffusion
      {
        std::vector<double> weights;
        double w_sum = 0;
        for (auto nid : inside_neighbours) {
          HERMES_ASSERT(nid < mesh->cell_centers.size());
          auto _w = w(mesh->cell_centers[cell_id], mesh->cell_centers[nid]);
          weights.emplace_back(_w);
          w_sum += _w;
        }
        HERMES_ASSERT(w_sum > 0);
        // normalize
        for (auto &weight : weights)
          weight /= w_sum;
        double l_sum = 0;
        for (size_t i = 0; i < inside_neighbours.size(); ++i)
          l_sum += weights[i] *
                   (field_values[cell_id] - field_values[inside_neighbours[i]]);
        s_field[cell_id] = l_sum;
      }
    }
  }
  return std::move(s_field);
}

} // namespace psa_anim
