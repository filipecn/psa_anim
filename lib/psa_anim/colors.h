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
///\file colors.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-12-13
///
///\brief

#ifndef PSA_ANIM_PSA_ANIM_COLORS_H
#define PSA_ANIM_PSA_ANIM_COLORS_H

#include <hermes/numeric/interpolation.h>
#include <hermes/numeric/numeric.h>
#include <vector>

namespace psa_anim {

struct Color {
  /// \brief Constructs color from rgba [0,255] values
  /// \param r
  /// \param g
  /// \param b
  /// \param a
  /// \return
  static Color fromU32(u32 r, u32 g, u32 b, u32 a = 255) {
    return {r / 255.f, g / 255.f, b / 255.f, a / 255.f};
  }
  /// \brief Extracts RGB components from unsigned integer
  /// \param color
  /// \return
  static Color rbgFromU32(u32 color) {
    return fromU32((color >> 16) & 0xFF, (color >> 8) & 0xFF, color & 0xFF);
  }
  /// \brief Default constructor (color black)
  Color() {
    r = g = b = 0.f;
    a = 1.f;
  }
  /// \brief 3 [0,1] float component constructor
  /// \param v
  explicit Color(const hermes::vec3 &v) : r(v.x), g(v.y), b(v.z), a(1.f) {}
  /// \brief rgba [0,1] components constructor
  /// \param _r
  /// \param _g
  /// \param _b
  /// \param _a
  Color(float _r, float _g, float _b, float _a = 1.f)
      : r(_r), g(_g), b(_b), a(_a) {}
  float r, g, b, a;
};

inline Color mix(float t, const Color &a, const Color &b) {
  return {hermes::interpolation::lerp(t, a.r, b.r),
          hermes::interpolation::lerp(t, a.g, b.g),
          hermes::interpolation::lerp(t, a.b, b.b)};
}

// *********************************************************************************************************************
//                                                                                                       ColorPalette
// *********************************************************************************************************************
/// \brief Color Palette representation
class ColorPalette {
public:
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  /// \brief Empty color palette constructor
  ColorPalette();
  /// \brief u32 data constructor
  /// \param c
  /// \param n
  explicit ColorPalette(const u32 *c, size_t n);
  /// \brief f32 data constructor
  /// \param c
  /// \param rgb_count
  explicit ColorPalette(const f32 *c, size_t rgb_count);
  /// \brief integer list data constructor
  /// \param c
  ColorPalette(std::initializer_list<int> c);
  /// \brief float list data constructor
  /// \param c
  ColorPalette(std::initializer_list<double> c);
  // *******************************************************************************************************************
  //                                                                                                        OPERATORS
  // *******************************************************************************************************************
  /// \brief Get color from parametric coordinate
  /// \param t
  /// \param alpha
  /// \return
  inline Color operator()(float t, float alpha = -1) const {
    t = hermes::Numbers::clamp(t, 0.f, 1.f);
    float ind = hermes::interpolation::lerp(
        t, 0.f, static_cast<float>(colors.size() - 1));
    float r = std::abs(hermes::Numbers::fract(ind));
    Color c;

    auto upper = hermes::Numbers::ceil2Int(ind);
    auto lower = hermes::Numbers::floor2Int(ind);

    if (upper >= static_cast<int>(colors.size()))
      c = colors[colors.size() - 1];
    else if (lower < 0)
      c = colors[0];
    else if (lower == upper)
      c = colors[lower];
    else
      c = mix(r, colors[lower], colors[upper]);
    if (alpha >= 0.0)
      c.a = alpha;
    else
      c.a = a;
    return c;
  }
  // *******************************************************************************************************************
  //                                                                                                    PUBLIC FIELDS
  // *******************************************************************************************************************
  float a{1.f};              //!< default alpha value
  std::vector<Color> colors; //!< raw color data
};

class TrigonometricColorPalette {
public:
  TrigonometricColorPalette(const hermes::vec3 &a, const hermes::vec3 &b,
                            const hermes::vec3 &c, const hermes::vec3 &d)
      : a{a}, b{b}, c{c}, d{d} {}
  inline Color operator()(float t, float alpha = 0) const {
    return Color(
        a[0] + b[0] * std::cos(hermes::Constants::two_pi * (c[0] * t + d[0])),
        a[1] + b[1] * std::cos(hermes::Constants::two_pi * (c[1] * t + d[1])),
        a[2] + b[2] * std::cos(hermes::Constants::two_pi * (c[1] * t + d[2])),
        alpha);
  }

  hermes::vec3 a;
  hermes::vec3 b;
  hermes::vec3 c;
  hermes::vec3 d;
};

struct ColorPalettes {
  // *******************************************************************************************************************
  //                                                                                                   STATIC METHODS
  // *******************************************************************************************************************
  /// \brief Matlab Heat Map color map
  /// \return
  static ColorPalette MatlabHeatMap();
  /// \brief Batlow color map
  /// \return
  static ColorPalette Batlow();

  static TrigonometricColorPalette Inigo6();
};

} // namespace psa_anim

#endif // PSA_ANIM_PSA_ANIM_COLORS_H
