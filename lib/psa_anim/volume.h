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
///\file volume.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2024-05-26
///
///\brief

#ifndef PSA_ANIM_VOLUME_H
#define PSA_ANIM_VOLUME_H

#include <hermes/geometry/bbox.h>
#include <hermes/logging/logging.h>
#include <mutex>
#include <openvdb/math/Transform.h>
#include <rasterizer.h>
#include <shared_mutex>

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <thread>
#include <unordered_map>

namespace psa_anim {

template <typename R> class Volume {
private:
  struct RasterizationGrid {
    std::size_t frame;
    std::vector<size_t> interpolation_regions;
    openvdb::FloatGrid::Ptr density;
  };

public:
  Volume(const PolyMesh &mesh, f32 voxel_size, f32 region_size,
         const hermes::vec3 &region_offset = {})
      : rasterizer_{R(mesh, voxel_size, region_size, region_offset)},
        transform_{
            openvdb::math::Transform::createLinearTransform(voxel_size)} {}

  R &rasterizer() { return rasterizer_; }

  const openvdb::math::Transform::Ptr transform() const { return transform_; }

  std::vector<size_t>
  fill(size_t frame, f32 min_value, const std::vector<double> &values,
       const std::vector<hermes::vec3> &grad_values,
       const std::function<bool(size_t, const hermes::point3 &, double &)>
           &voxel_callback = nullptr) {
    // filter cells
    std::vector<size_t> cell_ids;
    for (size_t cell_id = 0; cell_id < values.size(); ++cell_id)
      if (values[cell_id] >= min_value) {
        cell_ids.emplace_back(cell_id);
      }

    HERMES_LOG("filling {} cells", cell_ids.size());

    // prepare rasterizer
    // we will create/overwrite a vdb grid for each raster region
    auto region_ids = rasterizer_.setup(cell_ids);
    // prepare vdb grids
    rasterized_grids_.clear();
    size_t n_threads = std::thread::hardware_concurrency();
    for (size_t i = 0; i <= n_threads; ++i) {
      rasterized_grids_.emplace_back();
      rasterized_grids_.back().frame = frame;
      rasterized_grids_.back().density = openvdb::FloatGrid::create();
      rasterized_grids_.back().density->setTransform(transform_);
    }

    // fill regions by calling the rasterizer
    std::shared_mutex grids_mutex;
    size_t next_grid_index = 0;
    auto rasterized_regions = rasterizer_.raster(
        cell_ids, values, grad_values, [&](const VoxelCallbackData &data) {
          size_t grid_index = 0;
          {
            std::unique_lock<std::shared_mutex> lock(grids_mutex);
            grid_index = next_grid_index;
            next_grid_index++;
          }
          if (grid_index >= rasterized_grids_.size()) {
            HERMES_LOG_ERROR("{} >= {} -> {}", grid_index,
                             rasterized_grids_.size(), n_threads);
          }
          HERMES_ASSERT(grid_index < rasterized_grids_.size());
          openvdb::FloatGrid::Accessor density =
              rasterized_grids_[grid_index].density->getAccessor();
          size_t count = 0;
          for (const auto &item : data.regions) {
            HERMES_ASSERT(item.second.points.size() ==
                          item.second.values.size());
            HERMES_ASSERT(item.second.points.size() ==
                          item.second.indices.size());
            for (size_t i = 0; i < item.second.points.size(); ++i) {
              auto v = item.second.values[i];
              if (voxel_callback) {
                if (voxel_callback(item.first, item.second.points[i], v))
                  break;
              }
              count++;
              density.setValue({item.second.indices[i].i,
                                item.second.indices[i].j,
                                item.second.indices[i].k},
                               v);
            }
          }
        });
    HERMES_LOG_VARIABLE(next_grid_index);
    return rasterized_regions;
  }

  struct Local {
    static inline void write(const float &a, const float &b, float &result) {
      result = b;
    }
    static inline void max(openvdb::CombineArgs<float> &args) {
      if (args.b() > args.a()) {
        // Transfer the B value and its active state.
        args.setResult(args.b());
        args.setResultIsActive(args.bIsActive());
      } else {
        // Preserve the A value and its active state.
        args.setResult(args.a());
        args.setResultIsActive(args.aIsActive());
      }
    }
  };

  openvdb::FloatGrid::Ptr output(size_t frame) {
    HERMES_LOG("Merging {} rasterized grids in {}", rasterized_grids_.size(),
               std::this_thread::get_id());
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
    grid->setTransform(transform_);
    size_t i = 0;
    for (auto &rasterized_grid : rasterized_grids_) {
      grid->tree().combineExtended(rasterized_grid.density->tree(), Local::max);
    }
    return std::move(grid);
  }

  std::vector<size_t> interpolationRegions(size_t frame) {
    std::vector<size_t> region_ids;
    for (const auto &item : rasterized_grids_) {
      if (item.second.frame == frame)
        region_ids.emplace_back(item.first);
    }
    return region_ids;
  }

  hermes::bbox3 interpolationBox(size_t region_id) {
    return rasterizer_.interpolationGrid().regionBox(region_id);
  }

private:
  std::vector<RasterizationGrid> rasterized_grids_;
  R rasterizer_;
  openvdb::math::Transform::Ptr transform_;
};

} // namespace psa_anim

#endif // PSA_ANIM_VOLUME_H
