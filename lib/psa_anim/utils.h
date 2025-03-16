// Copyright (c) 2024, FilipeCN.
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
///\file utils.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2024-05-26
///
///\brief

#ifndef PSA_ANIM_UTILS_H
#define PSA_ANIM_UTILS_H

#include "hermes/common/file_system.h"
#include "hermes/common/result.h"
#include "hermes/geometry/transform.h"
#include "hermes/numeric/numeric.h"
#include <colors.h>
#include <openfoam_poly_mesh.h>

#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <queue>
#include <sys/resource.h>
#include <thread>

inline void RAM() {
  rusage ram;
  if (!getrusage(RUSAGE_SELF, &ram)) {
    HERMES_LOG("RAM USAGE: {} GB", ram.ru_maxrss / 1000000.);
  }
}

namespace psa_anim {

class ThreadPool {
public:
  ThreadPool(size_t num_threads = std::thread::hardware_concurrency()) {
    for (size_t i = 0; i < num_threads; ++i)
      threads_.emplace_back([this] {
        while (true) {
          std::function<void()> task;
          {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            cv_.wait(lock, [this] { return !tasks_.empty() || stop_; });

            if (stop_ && tasks_.empty())
              return;
            // Get the next task from the queue
            task = std::move(tasks_.front());
            tasks_.pop();
          }
          task();
        }
      });
  }
  ~ThreadPool() {
    {
      // Lock the queue to update the stop flag safely
      std::unique_lock<std::mutex> lock(queue_mutex_);
      stop_ = true;
    }

    // Notify all threads
    cv_.notify_all();

    // Joining all worker threads to ensure they have
    // completed their tasks
    for (auto &thread : threads_) {
      thread.join();
    }
  }

  // Enqueue task for execution by the thread pool
  void enqueue(std::function<void()> task) {
    {
      std::unique_lock<std::mutex> lock(queue_mutex_);
      tasks_.emplace(std::move(task));
    }
    cv_.notify_one();
  }

private:
  std::vector<std::thread> threads_;
  std::queue<std::function<void()>> tasks_;
  std::mutex queue_mutex_;
  std::condition_variable cv_;
  bool stop_ = false;
};

inline void for_indices(size_t n, size_t m,
                        const std::function<void(size_t, size_t)> &f) {
  if (n == 0)
    return;
  if (m == 0 || n <= m) {
    f(0, n - 1);
    return;
  }
  ThreadPool pool;
  for (size_t start = 0; start < n; start += m)
    pool.enqueue([start, n, m, f] {
      size_t end = std::min(start + m - 1, n - 1);
      f(start, end);
    });
}

inline hermes::Result<PolyMesh>
create_openfoam_grid_mesh(const hermes::size3 &res,
                          const hermes::vec3 &cell_size, float angle) {
  std::vector<hermes::point3> vertices;
  std::vector<std::vector<size_t>> faces;
  std::vector<std::vector<size_t>> cells;
  std::vector<PolyMesh::Patch> patches;

  // vertex flat index -> 3 min face indices (XY XZ YZ)
  // -1 indicates the vertex does not point a face in that direction
  std::vector<std::vector<int>> vertex_min_faces_map;

  hermes::vec3 domain_size{res.width * cell_size.x, res.depth * cell_size.y,
                           res.height * cell_size.z};

  hermes::range3 cell_range(res);
  hermes::range3 vertex_range(res + hermes::size3{1, 1, 1});
  for (auto ijk : vertex_range) {
    // add vertex

    float tan_angle_rcp =
        1 / std::tan(hermes::Trigonometry::degrees2radians(angle));
    hermes::point3 p = hermes::point3(ijk) * cell_size;
    hermes::vec3 shear(p.z * tan_angle_rcp, 0, 0);
    hermes::vec3 offset((domain_size.x + domain_size.z * tan_angle_rcp) * 0.5f,
                        domain_size.y * 0.5f, domain_size.z * 0.5f);

    vertices.emplace_back(p + shear - offset);
    //     _______
    //    /|_____/
    //   | |____|__
    //   | /    | /
    //   |/     |/
    //   --------
    //
    // add min faces
    std::vector<int> min_faces = {-1, -1, -1};
    // XY
    if (ijk.i < res.width && ijk.j < res.height) {
      min_faces[0] = faces.size();
      faces.push_back({vertex_range.flatIndex(ijk + hermes::index3(0, 0, 0)),
                       vertex_range.flatIndex(ijk + hermes::index3(1, 0, 0)),
                       vertex_range.flatIndex(ijk + hermes::index3(1, 1, 0)),
                       vertex_range.flatIndex(ijk + hermes::index3(0, 1, 0))});
    }
    // XZ
    if (ijk.i < res.width && ijk.k < res.depth) {
      min_faces[1] = faces.size();
      faces.push_back({vertex_range.flatIndex(ijk + hermes::index3(0, 0, 0)),
                       vertex_range.flatIndex(ijk + hermes::index3(1, 0, 0)),
                       vertex_range.flatIndex(ijk + hermes::index3(1, 0, 1)),
                       vertex_range.flatIndex(ijk + hermes::index3(0, 0, 1))});
    }
    // YZ
    if (ijk.j < res.height && ijk.k < res.depth) {
      min_faces[2] = faces.size();
      faces.push_back({vertex_range.flatIndex(ijk + hermes::index3(0, 0, 0)),
                       vertex_range.flatIndex(ijk + hermes::index3(0, 1, 0)),
                       vertex_range.flatIndex(ijk + hermes::index3(0, 1, 1)),
                       vertex_range.flatIndex(ijk + hermes::index3(0, 0, 1))});
    }

    vertex_min_faces_map.emplace_back(min_faces);
  }

  for (auto ijk : cell_range) {
    auto min_vertex_ijk = vertex_range.flatIndex(ijk);

    cells.push_back({
        vertex_min_faces_map[vertex_range.flatIndex(
            ijk + hermes::index3(0, 0, 0))][0],
        vertex_min_faces_map[vertex_range.flatIndex(
            ijk + hermes::index3(0, 0, 0))][1],
        vertex_min_faces_map[vertex_range.flatIndex(
            ijk + hermes::index3(0, 0, 0))][2],
        vertex_min_faces_map[vertex_range.flatIndex(
            ijk + hermes::index3(1, 0, 0))][2], // YZ
        vertex_min_faces_map[vertex_range.flatIndex(
            ijk + hermes::index3(0, 1, 0))][1], // XZ
        vertex_min_faces_map[vertex_range.flatIndex(
            ijk + hermes::index3(0, 0, 1))][0] // XY
    });
  }

  return PolyMesh::from(vertices, faces, cells, patches);
}

// procedural density tests

class Func {
public:
  virtual double f(const hermes::point3 &p) const = 0;
  virtual hermes::vec3 g(const hermes::point3 &p) const = 0;
  virtual std::vector<double> F(const std::vector<hermes::point3> &points) {
    std::vector<double> values;
    for (const auto &p : points)
      values.emplace_back(f(p));
    return values;
  }
  virtual std::vector<hermes::vec3>
  G(const std::vector<hermes::point3> &points) {
    std::vector<hermes::vec3> values;
    for (const auto &p : points)
      values.emplace_back(g(p));
    return values;
  }
};

class Gaussian : public Func {
public:
  virtual double f(const hermes::point3 &p) const override {
    auto x2 = hermes::distance2(p, center);
    return a * std::exp(-x2 / (2 * c * c));
  }
  virtual hermes::vec3 g(const hermes::point3 &p) const override {
    auto d = p - center;
    return d * (-(a * f(p)) / (c * c));
  }
  double a{1.0};
  double c{5.0};
  hermes::point3 center;
};
/*
class Gaussian : public Func {
public:
  virtual double f(const hermes::point3 &p) const override {
    static double k = 1. / std::sqrtf(hermes::Constants::two_pi);
    auto x2 = hermes::distance2(p, mean);
    return (1. / standard_deviation) * k *
           std::exp(-x2 / (2 * standard_deviation * standard_deviation));
  }
  virtual hermes::vec3 g(const hermes::point3 &p) const override {
    static double k = 1. / std::sqrtf(hermes::Constants::two_pi);
    auto d = mean - p;
    return d * k *
           (1. /
            (standard_deviation * standard_deviation * standard_deviation)) *
           f(p);
  }
  double standard_deviation{1.0};
  hermes::point3 mean;
};
*/

class Poly : public Func {
public:
  virtual double f(const hermes::point3 &p) const override {
    return p.x * p.x * p.x * p.x + 2 * p.y * p.y + p.z;
  }
  virtual hermes::vec3 g(const hermes::point3 &p) const override {
    return {4 * p.x * p.x * p.x, 4 * p.y, 1};
  }
};

inline void bbox2obj(const std::vector<hermes::bbox3> &boxes,
                     const hermes::Path &file) {
  hermes::Str s;
  // write vertices
  for (size_t i = 0; i < boxes.size(); ++i) {
    s.appendLine("o region.", i);
    const auto &box = boxes[i];
    for (size_t c = 0; c < 8; ++c) {
      auto v = box.corner(c);
      s.appendLine("v ", v.x, " ", v.y, " ", v.z);
    }
    size_t offset = i * 8 + 1;
    s.appendLine("l ", 0 + offset, " ", 1 + offset, " ", 3 + offset, " ",
                 2 + offset, " ", 0 + offset);
    s.appendLine("l ", 4 + offset, " ", 5 + offset, " ", 7 + offset, " ",
                 6 + offset, " ", 4 + offset);
    for (int j = 0; j < 4; ++j)
      s.appendLine("l ", j + offset, " ", j + 4 + offset);
  }
  hermes::FileSystem::writeFile(file, s.str());
}

FaceMesh polyMesh2faceMesh(const PolyMesh &mesh,
                           const std::vector<size_t> &cell_ids) {
  FaceMesh face_mesh;
  std::set<std::size_t> faces;
  for (auto cell_id : cell_ids) {
    auto cell_faces = mesh.cells[cell_id];
    for (auto cell_face : cell_faces)
      faces.insert(cell_face);
  }

  std::unordered_map<size_t, size_t> vertex_index;
  for (auto face_id : faces) {
    auto face_vertices = mesh.faces[face_id];
    std::vector<std::size_t> face;
    for (auto vertex_id : face_vertices) {
      if (!vertex_index.count(vertex_id)) {
        vertex_index[vertex_id] = face_mesh.vertices.size();
        face_mesh.vertices.emplace_back(mesh.vertices[vertex_id]);
      }
      face.emplace_back(vertex_index[vertex_id]);
    }
    face_mesh.faces.emplace_back(face);
  }

  return face_mesh;
}

FaceMesh polyMesh2faceMesh(const PolyMesh &mesh) {
  std::unordered_map<size_t, size_t> vertex_map;
  FaceMesh face_mesh;
  for (const auto &face : mesh.faces) {
    std::vector<size_t> face_vertices;
    for (auto v : face) {
      if (!vertex_map.count(v)) {
        vertex_map[v] = face_mesh.vertices.size();
        face_mesh.vertices.emplace_back(mesh.vertices[v]);
      }
      face_vertices.emplace_back(vertex_map[v]);
    }
    face_mesh.faces.emplace_back(face_vertices);
  }
  return face_mesh;
}

inline void bbox2ply(const std::vector<hermes::bbox3> &boxes,
                     const hermes::Path &file,
                     const std::vector<double> &values) {

  auto a_palette = psa_anim::ColorPalettes::Batlow();
  float max_a = 0;
  float min_a = 0;
  for (auto a : values) {
    max_a = std::max(max_a, (float)a);
    min_a = std::min(min_a, (float)a);
  }

  hermes::Str s;

  s.appendLine("ply");
  s.appendLine("format ascii 1.0");
  s.appendLine("element vertex ", 8 * boxes.size());
  s.appendLine("property float x");
  s.appendLine("property float y");
  s.appendLine("property float z");
  s.appendLine("property uchar red");
  s.appendLine("property uchar green");
  s.appendLine("property uchar blue");
  s.appendLine("element face ", 6 * boxes.size());
  s.appendLine("property list uchar int vertex_index");
  s.appendLine("end_header");

  for (size_t i = 0; i < boxes.size(); ++i) {
    // vertices
    auto box = boxes[i];
    hermes::vec3 center(box.center().x, box.center().y, box.center().z);
    f32 scale = 0.98;
    hermes::Transform t = hermes::Transform::translate(center) *
                          hermes::Transform::scale(scale, scale, scale) *
                          hermes::Transform::translate(-center);
    // box = t(box);
    auto color = a_palette((values[i] - min_a) / (max_a - min_a));
    for (size_t c = 0; c < 8; ++c) {
      auto v = box.corner(c);
      s.append(v.x, " ", v.y, " ", v.z);
      s.appendLine(" ", (int)(color.r * 255.f), " ", (int)(color.g * 255.f),
                   " ", (int)(color.b * 255.f));
    }
  }
  for (size_t i = 0; i < boxes.size(); ++i) {
    int faces[6][4] = {{0, 2, 3, 1}, {1, 3, 7, 5}, {4, 5, 7, 6},
                       {4, 6, 2, 0}, {2, 6, 7, 3}, {4, 0, 1, 5}};
    for (size_t j = 0; j < 6; ++j) {
      s.append("4");
      for (int k = 0; k < 4; ++k)
        s.append(" ", faces[j][k] + i * 8);
      s.appendLine();
    }
  }

  hermes::FileSystem::writeFile(file, s.str());
}

inline void points2ply(const std::vector<hermes::point3> &points,
                       const hermes::Path &file,
                       const std::vector<double> &values,
                       bool colormap = true) {

  auto a_palette = psa_anim::ColorPalettes::Batlow();
  float max_a = 0;
  float min_a = 0;
  for (auto a : values) {
    max_a = std::max(max_a, (float)a);
    min_a = std::min(min_a, (float)a);
  }

  hermes::Str s;

  s.appendLine("ply");
  s.appendLine("format ascii 1.0");
  s.appendLine("element vertex ", points.size());
  s.appendLine("property float x");
  s.appendLine("property float y");
  s.appendLine("property float z");
  s.appendLine("property uchar red");
  s.appendLine("property uchar green");
  s.appendLine("property uchar blue");
  s.appendLine("property list uchar int vertex_index");
  s.appendLine("end_header");

  for (size_t i = 0; i < points.size(); ++i) {
    const auto &v = points[i];
    // vertices
    s.append(v.x, " ", v.y, " ", v.z);
    if (colormap) {
      auto color = a_palette((values[i] - min_a) / (max_a - min_a));
      s.appendLine(" ", (int)(color.r * 255.f), " ", (int)(color.g * 255.f),
                   " ", (int)(color.b * 255.f));
    } else {
      auto t = (values[i] - min_a) / (max_a - min_a);
      s.appendLine(" ", (int)(t * 255.f), " ", (int)(t * 255.f), " ",
                   (int)(t * 255.f));
    }
  }

  hermes::FileSystem::writeFile(file, s.str());
}

inline void
stencils2ply(const std::vector<hermes::point3> &centers,
             const std::unordered_map<size_t, std::vector<size_t>> &stencils,
             const hermes::Path &file) {
  hermes::Str s;

  size_t vertex_count = 0;
  for (const auto &item : stencils)
    vertex_count += item.second.size();

  s.appendLine("ply");
  s.appendLine("format ascii 1.0");
  s.appendLine("element vertex ", vertex_count);
  s.appendLine("property float x");
  s.appendLine("property float y");
  s.appendLine("property float z");
  s.appendLine("property uchar red");
  s.appendLine("property uchar green");
  s.appendLine("property uchar blue");
  s.appendLine("property int region");
  s.appendLine("property int cell");
  s.appendLine("property list uchar int vertex_index");
  s.appendLine("end_header");

  auto palette = ColorPalettes::Inigo6();
  size_t region = 0;
  for (const auto &item : stencils) {
    auto color = palette(region++ / (float)stencils.size());
    for (const auto &cell_id : item.second) {
      const auto &v = centers[cell_id];
      s.append(v.x, " ", v.y, " ", v.z);
      s.append(" ", (int)(color.r * 255.f), " ", (int)(color.g * 255.f), " ",
               (int)(color.b * 255.f));
      s.append(" ", item.first);
      s.appendLine(" ", cell_id);
    }
  }
  hermes::FileSystem::writeFile(file, s.str());
}

} // namespace psa_anim

#endif
