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
///\file partition_grid.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2024-05-26
///
///\brief

#ifndef PSA_ANIM_PARITION_GRID_H
#define PSA_ANIM_PARITION_GRID_H

#include <hermes/logging/logging.h>
#include <openfoam_poly_mesh.h>
#include <utils.h>

#include <shared_mutex>
#include <unordered_map>
#include <vector>

namespace psa_anim {

class PartitionGrid {
public:
  PartitionGrid(const PolyMesh &mesh, f32 cell_size,
                const hermes::vec3 &offset = {})
      : mesh_{mesh} {
    f32 rcp = 1 / cell_size;
    toGrid_ = hermes::Transform::scale(rcp, rcp, rcp) *
              hermes::Transform::translate(-offset);
    toWorld_ = hermes::Transform::translate(offset) *
               hermes::Transform::scale(cell_size, cell_size, cell_size);
    // init cell map
    cell_partition_map_.resize(mesh_.cells.size());
    topology_ = psa_anim::computeFaceTopology(mesh.faces, mesh.cells);
  }

  const PolyMesh &mesh() const { return mesh_; }

  hermes::bbox3 regionBox(size_t key) const {
    if (key >= grid_positions_.size()) {
      HERMES_LOG_ERROR("key {} not found", key);
      return {};
    }
    auto lower = grid_positions_[key];
    auto upper = lower + hermes::index3(1, 1, 1);
    return {toWorld_(hermes::point3(lower.i, lower.j, lower.k)),
            toWorld_(hermes::point3(upper.i, upper.j, upper.k))};
  }

  size_t cellPartition(size_t cell_id) const {
    return cell_partition_map_[cell_id];
  }

  const std::vector<size_t> &regionCells(size_t region_id) const {
    static std::vector<size_t> dummy;
    if (region_id >= partition_cells_.size())
      return dummy;
    HERMES_ASSERT(region_id < partition_cells_.size());
    return partition_cells_[region_id];
  }

  const std::vector<size_t> &cellPartitions() const {
    return cell_partition_map_;
  }

  FaceMesh regionCellsMesh(size_t region_id) const {
    auto cells = regionCells(region_id);
    return psa_anim::polyMesh2faceMesh(mesh_, cells);
  }

  // setup partition regions
  std::set<size_t> init(const std::vector<size_t> cell_ids) {
    // retrieve grid hashes
    std::set<size_t> all_hashes;
    std::set<size_t> new_hashes;
    std::unordered_map<size_t, size_t> new_seeds;
    // touch new grid regions
    for (auto cell_id : cell_ids) {
      auto wp = mesh_.cell_centers[cell_id];
      size_t key = hash(wp);
      auto region = regionBox(key);
      if (!region.contains(wp)) {
        HERMES_LOG_VARIABLES(wp, region, key);
      }
      HERMES_ASSERT(region.contains(wp));
      all_hashes.insert(key);
      if (key < partition_cells_.size() && !partition_cells_[key].empty())
        continue;
      if (partition_cells_.size() <= key)
        partition_cells_.resize(key + 1);
      // register new hash
      new_hashes.insert(key);
      new_seeds[key] = cell_id;
    }
    // compute unitialized regions
    ThreadPool pool;
    for (auto key : new_hashes) {
      pool.enqueue([key, &new_seeds, this] {
        auto region = regionBox(key);
        partition_cells_[key] = psa_anim::intersect(
            topology_, mesh_.cells, mesh_.cell_centers, new_seeds[key], region);
        for (auto cell_id : partition_cells_[key])
          cell_partition_map_[cell_id] = key;
      });
    }
    return std::move(all_hashes);
  }

private:
  bool hash(const hermes::point3 &wp, size_t &key) const {
    auto gp = toGrid_(wp);
    hermes::index3 ijk(std::floor(gp.x), std::floor(gp.y), std::floor(gp.z));
    std::shared_lock<std::shared_mutex> lock(grid_positions_m_);
    auto i_it = hash_map_.find(ijk.i);
    if (i_it == hash_map_.end())
      return false;
    auto j_it = i_it->second.find(ijk.j);
    if (j_it == i_it->second.end())
      return false;
    auto k_it = j_it->second.find(ijk.k);
    if (k_it == j_it->second.end())
      return false;
    key = k_it->second;
    return true;
  }
  size_t hash(const hermes::point3 &wp) {
    auto gp = toGrid_(wp);
    hermes::index3 ijk(std::floor(gp.x), std::floor(gp.y), std::floor(gp.z));

    std::unordered_map<
        int, std::unordered_map<int, std::unordered_map<int, size_t>>>::iterator
        i_it = hash_map_.end();
    {
      std::shared_lock<std::shared_mutex> lock(grid_positions_m_);
      i_it = hash_map_.find(ijk.i);
    }
    if (i_it == hash_map_.end()) {
      std::unique_lock<std::shared_mutex> lock(grid_positions_m_);
      auto &hy = hash_map_[ijk.i][ijk.j];
      if (!hy.count(ijk.k)) {
        hy[ijk.k] = grid_positions_.size();
        grid_positions_.emplace_back(ijk);
        return grid_positions_.size() - 1;
      }
    }

    std::unordered_map<int, std::unordered_map<int, size_t>>::iterator j_it;
    {
      std::shared_lock<std::shared_mutex> lock(grid_positions_m_);
      j_it = i_it->second.find(ijk.j);
    }
    if (j_it == i_it->second.end()) {
      std::unique_lock<std::shared_mutex> lock(grid_positions_m_);
      auto &hy = i_it->second[ijk.j];
      if (!hy.count(ijk.k)) {
        hy[ijk.k] = grid_positions_.size();
        grid_positions_.emplace_back(ijk);
        return grid_positions_.size() - 1;
      }
    }

    std::unordered_map<int, size_t>::iterator k_it;
    {
      std::shared_lock<std::shared_mutex> lock(grid_positions_m_);
      k_it = j_it->second.find(ijk.k);
    }
    if (k_it == j_it->second.end()) {
      std::unique_lock<std::shared_mutex> lock(grid_positions_m_);
      j_it->second[ijk.k] = grid_positions_.size();
      grid_positions_.emplace_back(ijk);
      return grid_positions_.size() - 1;
    }

    return k_it->second;
  }

  // mesh
  const PolyMesh &mesh_;
  std::unordered_map<size_t, std::pair<i64, i64>> topology_;

  // openfoam cell -> partition
  std::vector<size_t> cell_partition_map_;

  // grid
  hermes::Transform toGrid_;
  hermes::Transform toWorld_;
  std::unordered_map<int,
                     std::unordered_map<int, std::unordered_map<int, size_t>>>
      hash_map_;
  std::vector<hermes::index3> grid_positions_;
  std::vector<std::vector<size_t>> partition_cells_;

  // parallel
  mutable std::shared_mutex grid_positions_m_;
};
} // namespace psa_anim

#endif
