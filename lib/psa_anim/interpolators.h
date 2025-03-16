/// Copyright (c) 2024, FilipeCN.
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
///\file interpolation.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2024-05-26
///
///\brief

#ifndef PSA_ANIM_TOOLS_PSA_ANIMPY_INTERPOLATION_GRID_H
#define PSA_ANIM_TOOLS_PSA_ANIMPY_INTERPOLATION_GRID_H

#include "utils.h"
#include <hermes/geometry/bbox.h>
#include <hermes/geometry/transform.h>
#include <hermes/logging/logging.h>
#include <interpolation.h>
#include <mutex>
#include <partition_grid.h>
#include <queue>
#include <set>
#include <thread>
#include <unordered_map>

namespace psa_anim {

// Nearest-point interpolation
// All points inside a cell receive the cell's value
class CellInterpolator {
public:
  CellInterpolator() = default;
  ~CellInterpolator() = default;
  void prepare(const PartitionGrid &grid, const std::set<size_t> &regions) {}

  void update(const PartitionGrid &grid, const std::set<size_t> regions,
              const std::vector<double> &values,
              const std::vector<hermes::vec3> &grad_values) {
    values_ = values.data();
  }

  std::vector<double>
  eval(const PartitionGrid &grid, size_t region,
       const std::vector<hermes::point3> &points,
       const std::vector<std::pair<size_t, size_t>> &cells) {
    std::vector<double> values;
    for (size_t i = 0; i < cells.size(); ++i)
      for (size_t j = 0; j < cells[i].second; ++j)
        values.emplace_back(values_[cells[i].first]);
    return values;
  }

private:
  const double *values_{nullptr};
};

class NearestInterpolator {
public:
  NearestInterpolator() = default;
  ~NearestInterpolator() = default;
  void prepare(const PartitionGrid &grid, const std::set<size_t> &regions) {
    std::vector<size_t> new_regions;
    for (auto region : regions) {
      if (!stencils_.count(region)) {
        stencils_[region] = {};
        new_regions.emplace_back(region);
      }
    }
    for_indices(new_regions.size(), 100,
                [this, &new_regions, &grid](size_t start, size_t end) {
                  for (size_t i = start; i <= end; ++i) {
                    size_t region = new_regions[i];
                    HERMES_ASSERT(stencils_.count(region));
                    Stencil stencil;

                    auto cell_ids = grid.regionCells(region);
                    for (auto cell_id : cell_ids) {
                      stencil.centers.emplace_back(
                          grid.mesh().cell_centers[cell_id]);
                      stencil.ids.emplace_back(cell_id);
                      stencil.values.emplace_back(0);
                    }
                    stencils_[region] = stencil;
                  }
                });
  }

  void update(const PartitionGrid &grid, const std::set<size_t> regions,
              const std::vector<double> &values,
              const std::vector<hermes::vec3> &grad_values) {
    std::vector<size_t> region_ids;
    for (auto region : regions)
      region_ids.emplace_back(region);
    for_indices(regions.size(), 100,
                [&region_ids, &grid, &values, &grad_values, this](size_t start,
                                                                  size_t end) {
                  for (size_t i = start; i <= end; ++i) {
                    HERMES_ASSERT(stencils_.count(region_ids[i]));
                    size_t region = region_ids[i];
                    auto &stencil = stencils_[region];
                    for (size_t j = 0; j < stencil.centers.size(); ++j)
                      stencil.values[j] = values[stencil.ids[j]];
                  }
                });
  }

  std::vector<double>
  eval(const PartitionGrid &grid, size_t region,
       const std::vector<hermes::point3> &points,
       const std::vector<std::pair<size_t, size_t>> &cells) {
    auto it = stencils_.find(region);
    HERMES_ASSERT(it != stencils_.end());
    Shepard shepard;
    std::vector<double> res;
    for (auto p : points)
      res.emplace_back(shepard.eval(it->second.centers, it->second.values, p));
    return res;
  }

private:
  struct Stencil {
    std::vector<hermes::point3> centers;
    std::vector<double> values;
    std::vector<size_t> ids;

    double closestValue(const hermes::point3 &p) {
      size_t closest_id = 0;
      float closest_dist = hermes::distance2(p, centers[0]);
      for (size_t i = 1; i < centers.size(); ++i) {
        float cur_dist = distance2(p, centers[i]);
        if (cur_dist < closest_dist) {
          closest_dist = cur_dist;
          closest_id = i;
        }
      }
      return values[closest_id];
    }
  };

  // region id -> (cell centers, values)
  std::unordered_map<size_t, Stencil> stencils_;
};

class ShepardInterpolator {
public:
  ShepardInterpolator() = default;
  ~ShepardInterpolator() = default;
  void prepare(const PartitionGrid &grid, const std::set<size_t> &regions) {
    std::vector<size_t> new_regions;
    for (auto region : regions) {
      if (!stencils_.count(region)) {
        stencils_[region] = {};
        new_regions.emplace_back(region);
      }
    }
    for_indices(new_regions.size(), 100,
                [this, &new_regions, &grid](size_t start, size_t end) {
                  for (size_t i = start; i <= end; ++i) {
                    size_t region = new_regions[i];
                    HERMES_ASSERT(stencils_.count(region));
                    Stencil stencil;

                    auto cell_ids = grid.regionCells(region);
                    for (auto cell_id : cell_ids) {
                      stencil.centers.emplace_back(
                          grid.mesh().cell_centers[cell_id]);
                      stencil.ids.emplace_back(cell_id);
                      stencil.values.emplace_back(0);
                    }
                    stencils_[region] = stencil;
                  }
                });
  }

  void update(const PartitionGrid &grid, const std::set<size_t> regions,
              const std::vector<double> &values,
              const std::vector<hermes::vec3> &grad_values) {
    std::vector<size_t> region_ids;
    for (auto region : regions)
      region_ids.emplace_back(region);
    for_indices(regions.size(), 1,
                [&region_ids, &grid, &values, &grad_values, this](size_t start,
                                                                  size_t end) {
                  for (size_t i = start; i <= end; ++i) {
                    HERMES_ASSERT(stencils_.count(region_ids[i]));
                    size_t region = region_ids[i];
                    auto &stencil = stencils_[region];
                    for (size_t j = 0; j < stencil.centers.size(); ++j) {
                      stencil.values[j] = values[stencil.ids[j]];
                    }
                  }
                });
  }

  std::vector<double>
  eval(const PartitionGrid &grid, size_t region,
       const std::vector<hermes::point3> &points,
       const std::vector<std::pair<size_t, size_t>> &cells) {
    auto it = stencils_.find(region);
    HERMES_ASSERT(it != stencils_.end());
    Shepard shepard;
    std::vector<double> res;
    for (auto p : points) {
      res.emplace_back(shepard.eval(it->second.centers, it->second.values, p));
    }
    return res;
  }

private:
  struct Stencil {
    std::vector<hermes::point3> centers;
    std::vector<double> values;
    std::vector<size_t> ids;
  };
  // region id -> (cell centers, values)
  std::unordered_map<size_t, Stencil> stencils_;
};

class KnnShepardInterpolator {
public:
  KnnShepardInterpolator() = default;
  ~KnnShepardInterpolator() = default;
  void prepare(const PartitionGrid &grid, const std::set<size_t> &regions) {
    std::vector<size_t> all_cells;
    // - gather all cells  with their centers
    for (auto region : regions) {
      auto cells = grid.regionCells(region);
      for (auto cell : cells) {
        cells_[cell].value = 0;
        cells_[cell].center = grid.mesh().cell_centers[cell];
        all_cells.emplace_back(cell);
      }
    }
    // - compute all knn for all cells
    for_indices(all_cells.size(), 100, [&](size_t start, size_t end) {
      hermes::point3 p;
      auto comp = [&](size_t a, size_t b) -> bool {
        return hermes::distance(cells_[a].center, p) >
               hermes::distance(cells_[b].center, p);
      };
      for (size_t i = start; i <= end; ++i) {
        auto this_cell = all_cells[i];
        p = cells_[this_cell].center;
        std::priority_queue<size_t, std::vector<size_t>, decltype(comp)> q(
            comp);
        for (auto cell_id : all_cells)
          q.push(cell_id);
        for (auto j = 0; j < KNN && !q.empty(); ++j) {
          stencils_[this_cell].emplace_back(q.top());
          q.pop();
        }
      }
    });
  }

  void update(const PartitionGrid &grid, const std::set<size_t> regions,
              const std::vector<double> &values,
              const std::vector<hermes::vec3> &grad_values) {
    // save all values
    for (auto &item : cells_)
      item.second.value = values[item.first];
  }

  std::vector<double>
  eval(const PartitionGrid &grid, size_t region,
       const std::vector<hermes::point3> &points,
       const std::vector<std::pair<size_t, size_t>> &cells) {
    Shepard shepard;
    std::vector<double> res;
    for (size_t i = 0; i < cells.size(); ++i) {
      auto item = cells[i];
      std::vector<hermes::point3> centers;
      std::vector<double> f;
      for (auto cell : stencils_[item.second]) {
        centers.emplace_back(cells_[cell].center);
        f.emplace_back(cells_[cell].value);
      }
      res.emplace_back(shepard.eval(centers, f, points[i]));
    }

    return res;
  }

private:
  struct Cell {
    double value;
    hermes::point3 center;
  };
  // cell id -> k closest cells
  std::unordered_map<size_t, std::vector<size_t>> stencils_;
  // cell id -> value
  std::unordered_map<size_t, Cell> cells_;
  size_t KNN{10};
};

// the hrbf is based on the cell_id, each cell belongs to a
// pre-defined region (interpolation grid cell)
// here we create the interpolation system (if necessary) by
// providing the cells that compose the respective region
template <typename RBF = CubicRBF> class HRBFInterpolator {
public:
  HRBFInterpolator() = default;
  ~HRBFInterpolator() = default;

  void prepare(const PartitionGrid &grid, const std::set<size_t> &regions) {
    std::vector<size_t> new_regions;
    for (auto region : regions) {
      if (!interpolants_.count(region)) {
        interpolants_[region] = {};
        new_regions.emplace_back(region);
      }
    }
    HERMES_LOG("new regions {}", new_regions.size());
    // compute unitialized regions
    for_indices(new_regions.size(), 100,
                [this, &new_regions, &grid](size_t start, size_t end) {
                  for (size_t i = start; i <= end; ++i) {
                    size_t region = new_regions[i];
                    interpolants_[region].init(grid.mesh().cell_centers,
                                               grid.regionCells(region));
                  }
                });
  }

  // build A
  void update(const PartitionGrid &grid, const std::set<size_t> regions,
              const std::vector<double> &values,
              const std::vector<hermes::vec3> &grad_values) {
    std::vector<size_t> region_ids;
    for (auto region : regions)
      region_ids.emplace_back(region);
    for_indices(regions.size(), 100,
                [&region_ids, &grid, &values, &grad_values, this](size_t start,
                                                                  size_t end) {
                  for (size_t i = start; i <= end; ++i) {
                    HERMES_ASSERT(interpolants_.count(region_ids[i]));
                    size_t region = region_ids[i];
                    interpolants_[region].build(values, grad_values,
                                                grid.regionCells(region));
                  }
                });
  }

  std::vector<double>
  eval(const PartitionGrid &grid, size_t region,
       const std::vector<hermes::point3> &points,
       const std::vector<std::pair<size_t, size_t>> &cells = {}) {
    HERMES_ASSERT(interpolants_.count(region));
    return interpolants_[region].eval(grid.mesh().cell_centers,
                                      grid.regionCells(region), points);
  }

private:
  std::unordered_map<size_t, InterpolationSystem<RBF>> interpolants_;
};

template <typename RBF = CubicRBF> class RBFInterpolator {
public:
  RBFInterpolator() = default;
  ~RBFInterpolator() = default;

  void prepare(const PartitionGrid &grid, const std::set<size_t> &regions) {
    std::vector<size_t> new_regions;
    for (auto region : regions) {
      if (!interpolants_.count(region)) {
        interpolants_[region] = {};
        new_regions.emplace_back(region);
      }
    }
    // compute unitialized regions
    for_indices(new_regions.size(), 100,
                [this, &new_regions, &grid](size_t start, size_t end) {
                  for (size_t i = start; i <= end; ++i) {
                    size_t region = new_regions[i];
                    interpolants_[region].init(grid.mesh().cell_centers,
                                               grid.regionCells(region));
                  }
                });
  }

  // build A
  void update(const PartitionGrid &grid, const std::set<size_t> regions,
              const std::vector<double> &values,
              const std::vector<hermes::vec3> &grad_values) {
    std::vector<size_t> region_ids;
    for (auto region : regions)
      region_ids.emplace_back(region);
    for_indices(regions.size(), 100,
                [&region_ids, &grid, &values, &grad_values, this](size_t start,
                                                                  size_t end) {
                  for (size_t i = start; i <= end; ++i) {
                    HERMES_ASSERT(interpolants_.count(region_ids[i]));
                    size_t region = region_ids[i];
                    interpolants_[region].build(values,
                                                grid.regionCells(region));
                  }
                });
  }

  std::vector<double>
  eval(const PartitionGrid &grid, size_t region,
       const std::vector<hermes::point3> &points,
       const std::vector<std::pair<size_t, size_t>> &cells = {}) {
    HERMES_ASSERT(interpolants_.count(region));
    return interpolants_[region].eval(grid.mesh().cell_centers,
                                      grid.regionCells(region), points);
  }

private:
  std::unordered_map<size_t, RBFSystem<RBF>> interpolants_;
};

} // namespace psa_anim

#endif // PSA_ANIM_TOOLS_PSA_ANIMPY_INTERPOLATION_GRID_H
