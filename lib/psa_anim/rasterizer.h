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

#ifndef PSA_ANIM_TOOLS_PSA_ANIMPY_FIELD_SAMPLER_H
#define PSA_ANIM_TOOLS_PSA_ANIMPY_FIELD_SAMPLER_H

#include "utils.h"
#include <hermes/geometry/vector.h>
#include <hermes/logging/logging.h>
#include <hermes/numeric/numeric.h>
#include <mutex>
#include <openfoam_poly_mesh.h>
#include <partition_grid.h>
#include <thread>
#include <unordered_map>

namespace psa_anim {
struct VoxelCallbackData {
  struct VoxelRegionData {
    std::vector<hermes::point3> points;
    std::vector<hermes::index3> indices;
    std::vector<double> values;
  };
  // interpolation region -> data
  std::unordered_map<size_t, VoxelRegionData> regions;
};

template <typename Interpolator> class Rasterizer {
public:
  Rasterizer(const PolyMesh &mesh, f32 voxel_size, f32 region_size,
             const hermes::vec3 &region_offset = {})
      : grid_{mesh, region_size, region_offset}, voxel_size_(voxel_size) {}

  std::set<size_t> setup(const std::vector<size_t> &cell_ids) {
    auto regions = grid_.init(cell_ids);
    HERMES_PING;
    interpolator_.prepare(grid_, regions);
    return regions;
  }

  const Interpolator &interpolator() const { return interpolator_; }

  const PartitionGrid &interpolationGrid() const { return grid_; }

  void set2d(bool is_2d) { is_2d_ = is_2d; }

  void setOnlyCells(bool only_cells) { rasterize_only_cells_ = only_cells; }

  std::vector<size_t> raster(
      const std::vector<size_t> &cell_ids, const std::vector<double> &values,
      const std::vector<hermes::vec3> &grad_values,
      const std::function<void(const VoxelCallbackData &)> &voxels_callback) {
    if (cell_ids.empty())
      return {};
    // prepare cell sets (one per partition)
    // partition -> cells
    std::unordered_map<size_t, std::vector<size_t>> partition_cells;
    for (auto cell_id : cell_ids) {
      auto region_id = grid_.cellPartition(cell_id);
      partition_cells[region_id].emplace_back(cell_id);
    }
    // prepare interpolator
    std::set<size_t> region_ids;
    for (const auto &item : partition_cells)
      region_ids.insert(item.first);
    HERMES_LOG("updating interpolator with {} regions", region_ids.size());
    interpolator_.update(grid_, region_ids, values, grad_values);
    HERMES_LOG("...done");
    // auxiliary data
    std::vector<size_t> flat_region_ids;
    for (size_t id : region_ids)
      flat_region_ids.emplace_back(id);
    const auto &mesh = grid_.mesh();
    auto w2v = hermes::Transform::scale(1 / voxel_size_, 1 / voxel_size_,
                                        1 / voxel_size_);
    auto v2w = hermes::inverse(w2v);
    size_t n_threads = std::thread::hardware_concurrency();
    HERMES_LOG_WARNING("rasterizing {} regions with {} threads ({} per thread)",
                       flat_region_ids.size(), n_threads,
                       flat_region_ids.size() / n_threads);
    if (!rasterize_only_cells_) {
      for_indices(
          flat_region_ids.size(),
          std::ceil((double)flat_region_ids.size() / (double)n_threads),
          [&](size_t start, size_t end) {
            VoxelCallbackData callback_data;
            for (size_t i = start; i <= end; ++i) {
              auto partition_id = flat_region_ids[i];
              auto &partition_data = callback_data.regions[partition_id];
              auto region = grid_.regionBox(partition_id);

              auto v_region = w2v(region);
              hermes::range3 range(
                  hermes::index3(hermes::Numbers::floor2Int(v_region.lower.x),
                                 hermes::Numbers::floor2Int(v_region.lower.y),
                                 hermes::Numbers::floor2Int(v_region.lower.z)),
                  hermes::index3(hermes::Numbers::ceil2Int(v_region.upper.x),
                                 hermes::Numbers::ceil2Int(v_region.upper.y),
                                 hermes::Numbers::ceil2Int(v_region.upper.z)) +
                      hermes::index3(1, 1, 1));

              for (auto ijk : range) {
                auto wp = v2w(hermes::point3(ijk.i, ijk.j, ijk.k));
                if (region.contains(wp) &&
                    ((is_2d_ && ijk.j == 0) || !is_2d_)) {
                  partition_data.points.emplace_back(wp);
                  partition_data.indices.emplace_back(ijk);
                }
              }

              // interpolate
              if (!partition_data.points.empty())
                partition_data.values = interpolator_.eval(
                    grid_, partition_id, partition_data.points, {});
            }
            voxels_callback(callback_data);
          });
    } else {
      for_indices(
          flat_region_ids.size(),
          std::ceil((double)flat_region_ids.size() / (double)n_threads),
          [&](size_t start, size_t end) {
            // HERMES_LOG_VARIABLE(std::this_thread::get_id());
            VoxelCallbackData callback_data;
            for (size_t i = start; i <= end; ++i) {
              auto partition_id = flat_region_ids[i];
              auto &partition_data = callback_data.regions[partition_id];

              auto it = partition_cells.find(partition_id);
              HERMES_ASSERT(it != partition_cells.end());
              const auto &cell_ids = it->second;

              // cell id - voxel count
              std::vector<std::pair<size_t, size_t>> cell_voxels;

              for (size_t cell_id : cell_ids) {
                // Check voxels against polygon
                const auto &region = mesh.cell_bounds[cell_id];
                auto v_region = w2v(region);
                hermes::range3 range(
                    hermes::index3(
                        hermes::Numbers::floor2Int(v_region.lower.x),
                        hermes::Numbers::floor2Int(v_region.lower.y),
                        hermes::Numbers::floor2Int(v_region.lower.z)),
                    hermes::index3(
                        hermes::Numbers::ceil2Int(v_region.upper.x),
                        hermes::Numbers::ceil2Int(v_region.upper.y),
                        hermes::Numbers::ceil2Int(v_region.upper.z)) +
                        hermes::index3(1, 1, 1));

                size_t voxel_in_cell_count = 0;
                for (auto ijk : range) {
                  auto wp = v2w(hermes::point3(ijk.i, ijk.j, ijk.k));
                  if (mesh.cellContainsPoint(cell_id, wp)) {
                    partition_data.points.emplace_back(wp);
                    partition_data.indices.emplace_back(ijk);
                    voxel_in_cell_count++;
                  }
                }
                if (voxel_in_cell_count)
                  cell_voxels.emplace_back(
                      std::make_pair(cell_id, voxel_in_cell_count));
              }
              // interpolate
              partition_data.values = interpolator_.eval(
                  grid_, partition_id, partition_data.points, cell_voxels);
            }

            voxels_callback(callback_data);
          });
    }
    return flat_region_ids;
  }

private:
  bool rasterize_only_cells_{false};
  PartitionGrid grid_;
  f32 voxel_size_;
  Interpolator interpolator_;
  bool is_2d_{false};
};

} // namespace psa_anim

#endif
