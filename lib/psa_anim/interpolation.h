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
///\file boundary_marker.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-11-14
///
///\brief

#ifndef PSA_ANIM_TOOLS_PSA_ANIMPY_INTERPOLATION_H
#define PSA_ANIM_TOOLS_PSA_ANIMPY_INTERPOLATION_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <hermes/common/file_system.h>
#include <hermes/geometry/point.h>
#include <hermes/logging/logging.h>
#include <hermes/numeric/numeric.h>
#include <string>

#define POLY 4
// #define DEBUG

namespace psa_anim {

struct MultiQuadricRBF {
  static double phi(double r, double e) {
    return std::sqrt(1 + (e * r) * (e * r));
  }
  static double ddr(double r, double e) {
    return e * e * r / MultiQuadricRBF::phi(r, e);
  }
  static double d2dr2(double r, double e) {
    double _phi = MultiQuadricRBF::phi(r, e);
    return e * e / (_phi * _phi * _phi);
  }
  static double ddx(double dx, double r, double e) {
    return dx * e * e / MultiQuadricRBF::phi(r, e);
  }
  static double d2dx2(double dx, double r, double e) {
    double e2 = e * e;
    double phi2 = 1 + e2 * r * r;
    double _phi = std::sqrt(phi2);
    return e2 * (phi2 - dx * dx * e2) / (_phi * _phi * _phi);
  }
  static double d2dxy(double dx, double dy, double r, double e) {
    double _phi = MultiQuadricRBF::phi(r, e);
    double e2 = e * e;
    return -dx * dy * e2 * e2 / (_phi * _phi * _phi);
  }
};

struct CubicRBF {
  static double phi(double r, double e = 0) { return r * r * r; }
  static double ddr(double r, double e = 0) { return 3 * r * r; }
  static double d2dr2(double r, double e = 0) { return 6 * r; }
  static double ddx(double dx, double r, double e = 0) { return dx * 3 * r; }
  static double d2dx2(double dx, double r, double e = 0) {
    if (hermes::Check::is_zero(r))
      return 0;
    return 3 * (r + dx * dx / r);
  }
  static double d2dxy(double dx, double dy, double r, double e = 0) {
    if (hermes::Check::is_zero(r))
      return 0;
    return 3 * dx * dy / r;
  }
};

std::string toString(const Eigen::MatrixXd &M) {
  hermes::Str s;
  for (size_t i = 0; i < M.rows(); ++i) {
    for (size_t j = 0; j < M.cols(); ++j) {
      s.append(std::setw(10), M(i, j), " ");
    }
    s.appendLine();
  }
  return s.str();
}

void print(const Eigen::MatrixXd &M) { std::cout << toString(M) << std::endl; }

class MatrixS {
public:
  MatrixS(int n, int m) {
    data.resize(n);
    for (int i = 0; i < n; ++i)
      data[i].resize(m);
  }
  hermes::Str &operator()(int i, int j) { return data[i][j]; }
  const hermes::Str &operator()(int i, int j) const { return data[i][j]; }
  void print(const std::string &filename, const MatrixS &sign, int n) {
    hermes::Str s;
    for (size_t i = 0; i < data.size(); ++i) {
      if (i % n == 0)
        s.appendLine();
      for (size_t j = 0; j < data[i].size(); ++j) {
        if (j % n == 0)
          s.append("   ");
        s.append(std::setw(8), hermes::Str::concat(sign.data[i][j], data[i][j]),
                 "  ");
      }
      s.appendLine();
    }
    hermes::FileSystem::writeFile(filename, s.str());
  }
  std::vector<std::vector<hermes::Str>> data;
};

template <typename RBFKernel> class InterpolationSystem {
public:
  void init(const std::vector<hermes::point3> &centers,
            const std::vector<size_t> &ids, double e = 0) {
#define P(I) centers[ids[I]]
#define D(I, J) hermes::distance(P(I), P(J))

#define UPPER_diag BLOCK_I_OFFSET + i, BLOCK_J_OFFSET + i
#define UPPER_upper BLOCK_I_OFFSET + i, BLOCK_J_OFFSET + j
#define UPPER_lower BLOCK_I_OFFSET + j, BLOCK_J_OFFSET + i
#define LOWER_diag BLOCK_J_OFFSET + i, BLOCK_I_OFFSET + i
#define LOWER_upper BLOCK_J_OFFSET + i, BLOCK_I_OFFSET + j
#define LOWER_lower BLOCK_J_OFFSET + j, BLOCK_I_OFFSET + i
#define DIAG_diag BLOCK_I_OFFSET + i, BLOCK_I_OFFSET + i
#define DIAG_upper BLOCK_I_OFFSET + i, BLOCK_I_OFFSET + j
#define DIAG_lower BLOCK_I_OFFSET + j, BLOCK_I_OFFSET + i

    //                        p  px py pz
    //     N    N    N    N   1  1  1  1
    // N  PHI -DDX -DDY -DDZ  1  X  Y  Z      w     alpha
    // N  DDX -DXX -DXY -DXZ  0  1  0  0      wx    alpha_x
    // N  DDY -DXY -DYY -DYZ  0  0  1  0   x  wy =  alpha_y
    // N  DDZ -DXZ -DYZ -DZZ  0  0  0  1      wz    alpha_z
    // 1   1    0    0    0   0  0  0  0      g     0
    // 1   X    1    0    0   0  0  0  0      gx    0
    // 1   Y    0    1    0   0  0  0  0      gy    0
    // 1   Z    0    0    1   0  0  0  0      gz    0
    // setup matrix size
    size_t n = ids.size();
    size_t N = n * 4 + POLY;
    A = Eigen::MatrixXd::Zero(N, N);

#ifdef DEBUG
    MatrixS M(N, N); // entry name
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        M(i, j) = "0";
    MatrixS S(N, N); // entry signal
                     //
    char DIM[] = {'X', 'Y', 'Z'};
#endif
    // p
    if (POLY > 0) {
      for (size_t i = 0; i < n; ++i) {
        // 1
        int offset = n * 4;
        for (size_t dim = 0; dim <= 3; ++dim) {
          A(i + n * dim, offset + dim) = A(offset + dim, i + n * dim) = 1;
#ifdef DEBUG
          M(i + n * dim, offset + dim) = M(offset + dim, i + n * dim) = "1";
#endif
          // X  Y  Z
          if (dim > 0) {
            A(i, offset + dim) = A(offset + dim, i) = P(i)[dim - 1];
#ifdef DEBUG
            M(i, offset + dim) = M(offset + dim, i) = DIM[dim - 1];
#endif
          }
        }
      }
    }

    // PHI
    for (size_t i = 0; i < n; ++i) {
      A(i, i) = RBFKernel::phi(0, e);
#ifdef DEBUG
      M(i, i) = "PHI" + std::to_string((i + 1) * 10 + (i + 1));
#endif
      for (size_t j = i + 1; j < n; ++j) {
        A(i, j) = A(j, i) = RBFKernel::phi(D(i, j), e);
#ifdef DEBUG
        M(i, j) = M(j, i) = "PHI" + std::to_string((i + 1) * 10 + (j + 1));
#endif
      }
    }
    // DD?
    for (int dim = 0; dim < 3; ++dim) {
      int BLOCK_I_OFFSET = 0;
      int BLOCK_J_OFFSET = dim * n + n;
      for (int i = 0; i < n; ++i) {
        // DD? diag
        A(UPPER_diag) = -RBFKernel::ddx(0, 0, e);
        A(LOWER_diag) = -A(UPPER_diag);
#ifdef DEBUG
        S(UPPER_diag) = "-";
        S(LOWER_diag) = "+";
        M(UPPER_diag) = M(LOWER_diag) =
            hermes::Str::concat("DD", DIM[dim], "(", (i + 1) * 10 + i + 1, ")");
#endif

        for (int j = i + 1; j < n; ++j) {
          auto d = P(j) - P(i);

          // d/d dim (phi(i,j)) -> -DD?(I, J)
          A(UPPER_upper) = RBFKernel::ddx(d[dim], d.length(), e);
          A(UPPER_lower) = -A(UPPER_upper);

          // d/d dim (phi(i,j)) -> DD?(J, I) = -DD?(I, J)
          A(LOWER_upper) = A(UPPER_lower);
          A(LOWER_lower) = A(UPPER_upper);

#ifdef DEBUG
          S(UPPER_upper) = "-";
          S(UPPER_lower) = "+";
          S(LOWER_upper) = "+";
          S(LOWER_lower) = "-";
          M(UPPER_upper) = M(UPPER_lower) = M(LOWER_upper) = M(LOWER_lower) =
              hermes::Str::concat("DD", DIM[dim], "(", (i + 1) * 10 + j + 1,
                                  ")");
#endif
        }
      }
    }

    // DXX
    for (int dim = 0; dim < 3; ++dim) {
      int BLOCK_I_OFFSET = dim * n + n;
      for (int i = 0; i < n; ++i) {
        A(DIAG_diag) = -RBFKernel::d2dx2(0, 0, e);
#ifdef DEBUG
        M(DIAG_diag) = hermes::Str::concat("D", DIM[dim], DIM[dim], "(",
                                           (i + 1) * 10 + i + 1, ")");
        S(DIAG_diag) = "-";
#endif
        for (int j = i + 1; j < n; ++j) {
          // d/d dim (phi(i,j)) -> -DD?(I, J)
          auto d = P(j) - P(i);
          A(DIAG_upper) = A(DIAG_lower) =
              -RBFKernel::d2dx2(d[dim], d.length(), e);
#ifdef DEBUG
          M(DIAG_upper) = M(DIAG_lower) = hermes::Str::concat(
              "D", DIM[dim], DIM[dim], "(", (i + 1) * 10 + j + 1, ")");
          S(DIAG_upper) = S(DIAG_lower) = "-";
#endif
        }
      }
    }

    // DXY
    for (int dim_i = 0; dim_i < 3; ++dim_i) {
      int BLOCK_I_OFFSET = dim_i * n + n;
      for (int dim_j = dim_i + 1; dim_j < 3; ++dim_j) {
        int BLOCK_J_OFFSET = dim_j * n + n;
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            auto d = P(j) - P(i);

            A(UPPER_upper) = A(UPPER_lower) = A(LOWER_upper) = A(LOWER_lower) =
                -RBFKernel::d2dxy(d[dim_i], d[dim_j], d.length(), e);
#ifdef DEBUG
            M(UPPER_upper) = M(UPPER_lower) = M(LOWER_upper) = M(LOWER_lower) =
                hermes::Str::concat("D", DIM[dim_i], DIM[dim_j], "(",
                                    (i + 1) * 10 + j + 1, ")");
            S(UPPER_upper) = S(UPPER_lower) = S(LOWER_upper) = S(LOWER_lower) =
                "-";
#endif
          }
        }
      }
    }
#ifdef DEBUG
    M.print("dump", S, n);
    hermes::FileSystem::writeFile("A", toString(A));
#endif
    A = A.inverse();
    w.resize(N);
    b.resize(N);
#undef D
#undef P
  }

  void build(const std::vector<double> &values,
             const std::vector<hermes::vec3> &grad_values,
             const std::vector<size_t> &ids) {
#define F(I) values[ids[I]]
#define G(D, I) grad_values[ids[I]][D]

    size_t n = ids.size();
    size_t N = n * 4 + POLY;

    if (w.rows() != N || b.rows() != N) {
      HERMES_LOG_ERROR("bad system size {} != {} != {}", N, w.rows(), b.rows());
      return;
    }

    for (size_t i = 0; i < n; ++i)
      b[i] = F(i);
    for (size_t dim = 0; dim < 3; ++dim)
      for (size_t i = 0; i < n; ++i)
        b[n + dim * n + i] = G(dim, i);

    for (size_t p = 0; p < POLY; ++p)
      b[4 * n + p] = 0;

    w = A * b;

#undef G
#undef F
  }

  std::vector<double> eval(const std::vector<hermes::point3> &centers,
                           const std::vector<size_t> &ids,
                           const std::vector<hermes::point3> &points,
                           double e = 0) const {
#define C(I) centers[ids[I]]

    size_t n = ids.size();
    size_t N = n * 4 + POLY;

    if (w.rows() != N) {
      HERMES_LOG_ERROR("bad system size {} != {} != {}", N, w.rows(), b.rows());
      return {};
    }

    Eigen::VectorXd v;
    auto wt = w.transpose();
    v.resize(N);

    std::vector<double> r;

    bool first = true;
    for (const auto &p : points) {
      for (size_t i = 0; i < n; ++i)
        v[i] = RBFKernel::phi(hermes::distance(p, C(i)), e);
      for (size_t dim = 0; dim < 3; ++dim) {
        for (size_t i = 0; i < n; ++i) {
          auto d = p - C(i);
          v[n + dim * n + i] = -RBFKernel::ddx(d[dim], d.length(), e);
        }
        if (POLY == 4) {
          v[4 * n + 0] = 1;
          v[4 * n + dim + 1] = p[dim];
        }
      }

      r.emplace_back(wt * v);
    }

    return r;

#undef C
  }

  bool empty() const { return A.size() == 0; }

private:
  Eigen::MatrixXd A;
  Eigen::VectorXd w;
  Eigen::VectorXd b;
};

template <typename RBFKernel> class RBFSystem {
public:
  void init(const std::vector<hermes::point3> &centers,
            const std::vector<size_t> &ids, double e = 0) {
#define P(I) centers[ids[I]]
#define D(I, J) hermes::distance(P(I), P(J))

#define UPPER_diag BLOCK_I_OFFSET + i, BLOCK_J_OFFSET + i
#define UPPER_upper BLOCK_I_OFFSET + i, BLOCK_J_OFFSET + j
#define UPPER_lower BLOCK_I_OFFSET + j, BLOCK_J_OFFSET + i
#define LOWER_diag BLOCK_J_OFFSET + i, BLOCK_I_OFFSET + i
#define LOWER_upper BLOCK_J_OFFSET + i, BLOCK_I_OFFSET + j
#define LOWER_lower BLOCK_J_OFFSET + j, BLOCK_I_OFFSET + i
#define DIAG_diag BLOCK_I_OFFSET + i, BLOCK_I_OFFSET + i
#define DIAG_upper BLOCK_I_OFFSET + i, BLOCK_I_OFFSET + j
#define DIAG_lower BLOCK_I_OFFSET + j, BLOCK_I_OFFSET + i

    //     N
    // N  PHI    w   alpha

    // setup matrix size
    size_t n = ids.size();
    size_t N = n;
    A = Eigen::MatrixXd::Zero(N, N);

#ifdef DEBUG
    MatrixS M(N, N); // entry name
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        M(i, j) = "0";
    MatrixS S(N, N); // entry signal
                     //
    char DIM[] = {'X', 'Y', 'Z'};
#endif

    // PHI
    for (size_t i = 0; i < n; ++i) {
      A(i, i) = RBFKernel::phi(0, e);
#ifdef DEBUG
      M(i, i) = "PHI" + std::to_string((i + 1) * 10 + (i + 1));
#endif
      for (size_t j = i + 1; j < n; ++j) {
        A(i, j) = A(j, i) = RBFKernel::phi(D(i, j), e);
#ifdef DEBUG
        M(i, j) = M(j, i) = "PHI" + std::to_string((i + 1) * 10 + (j + 1));
#endif
      }
    }

#ifdef DEBUG
    M.print("dump", S, n);
    hermes::FileSystem::writeFile("A", toString(A));
#endif
    A = A.inverse();
    w.resize(N);
    b.resize(N);
#undef D
#undef P
  }

  void build(const std::vector<double> &values,
             const std::vector<size_t> &ids) {
#define F(I) values[ids[I]]
#define G(D, I) grad_values[ids[I]][D]

    size_t n = ids.size();
    size_t N = n;

    if (w.rows() != N || b.rows() != N) {
      HERMES_LOG_ERROR("bad system size {} != {} != {}", N, w.rows(), b.rows());
      return;
    }

    for (size_t i = 0; i < n; ++i)
      b[i] = F(i);

    w = A * b;

#undef G
#undef F
  }

  std::vector<double> eval(const std::vector<hermes::point3> &centers,
                           const std::vector<size_t> &ids,
                           const std::vector<hermes::point3> &points,
                           double e = 0) const {
#define C(I) centers[ids[I]]

    size_t n = ids.size();
    size_t N = n;

    if (w.rows() != N) {
      HERMES_LOG_ERROR("bad system size {} != {} != {}", N, w.rows(), b.rows());
      return {};
    }

    Eigen::VectorXd v;
    auto wt = w.transpose();
    v.resize(N);

    std::vector<double> r;

    bool first = true;
    for (const auto &p : points) {
      for (size_t i = 0; i < n; ++i)
        v[i] = RBFKernel::phi(hermes::distance(p, C(i)), e);

      r.emplace_back(wt * v);
    }

    return r;

#undef C
  }

  bool empty() const { return A.size() == 0; }

private:
  Eigen::MatrixXd A;
  Eigen::VectorXd w;
  Eigen::VectorXd b;
};

class Shepard {
public:
  double eval(const std::vector<hermes::point3> &centers,
              const std::vector<double> &f, const hermes::point3 &p) const {
    double den = 0.0;
    double num = 0.0;

    auto cmp = [&](int a, int b) -> bool {
      auto dista = hermes::distance(centers[a], p);
      auto distb = hermes::distance(centers[b], p);
      return dista > distb;
    };
    std::priority_queue<int, std::vector<int>, decltype(cmp)> q(cmp);
    for (int i = 0; i < centers.size(); ++i) {
      q.push(i);
    }
    int k = 0;
    while (!q.empty() && k < 10) {
      k++;
      auto i = q.top();
      q.pop();
      auto dist = hermes::distance(centers[i], p);
      double phi = 1.0 / (1e-6 + dist);
      num += f[i] * phi;
      den += phi;
    }
    // for (int i = 0; i < centers.size(); ++i) {
    //   auto dist = hermes::distance(centers[i], p);
    //   double phi = 1.0 / (1e-6 + dist);
    //   num += f[i] * phi;
    //   den += phi;
    // }
    HERMES_ASSERT(den != 0);
    return num / den;
  }
};

} // namespace psa_anim

#endif // PSA_ANIM_TOOLS_PSA_ANIMPY_INTERPOLATION_H
