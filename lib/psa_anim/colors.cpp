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
///\file colors.cpp
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-12-13
///
///\brief

#include <colors.h>

namespace psa_anim {

/// \brief Batlow color map data
/// \note Extracted from https://www.fabiocrameri.ch/colourmaps/
unsigned int batlow[] = {
    0x011959, 0x041E5A, 0x06215B, 0x08255B, 0x0A285C, 0x0B2D5D, 0x0C2F5D,
    0x0D335E, 0x0E375E, 0x0F395F, 0x103D5F, 0x103F60, 0x114260, 0x114460,
    0x124761, 0x134A61, 0x144C62, 0x154F62, 0x165062, 0x175362, 0x195662,
    0x1A5762, 0x1C5A62, 0x1E5C62, 0x205E61, 0x226061, 0x256260, 0x28645F,
    0x2B655E, 0x2F675C, 0x31695B, 0x356A59, 0x386C58, 0x3C6D56, 0x416F53,
    0x447052, 0x48714F, 0x4C724D, 0x50744B, 0x537549, 0x587646, 0x5D7844,
    0x607942, 0x657A3F, 0x687B3E, 0x6D7C3B, 0x717D39, 0x767F37, 0x7B8034,
    0x7F8133, 0x848331, 0x88842F, 0x8E852E, 0x93872C, 0x97882C, 0x9D892B,
    0xA18A2B, 0xA78B2C, 0xAB8C2D, 0xB18D2F, 0xB78E31, 0xBB8F33, 0xC09036,
    0xC49138, 0xCA923C, 0xCD923F, 0xD29343, 0xD89448, 0xDB954B, 0xE09651,
    0xE39754, 0xE7985A, 0xEA995E, 0xEE9B64, 0xF19D6B, 0xF39E70, 0xF6A077,
    0xF8A17B, 0xF9A382, 0xFAA587, 0xFCA78E, 0xFCA995, 0xFDAB9A, 0xFDADA0,
    0xFDAFA5, 0xFDB1AB, 0xFDB3B1, 0xFDB4B6, 0xFDB7BC, 0xFDB8C0, 0xFDBAC7,
    0xFDBCCB, 0xFCBED1, 0xFCC0D8, 0xFCC2DC, 0xFCC4E3, 0xFBC6E8, 0xFBC8EF,
    0xFBCAF3};

/// \brief Matlab Heat color map data
float matlab_heat[] = {
    0.64706f,  0.f,       0.14902f, 0.67059f,   0.015686f, 0.14902f,
    0.69412f,  0.035294f, 0.14902f, 0.71373f,   0.05098f,  0.14902f,
    0.73725f,  0.070588f, 0.14902f, 0.75686f,   0.086275f, 0.14902f,
    0.77647f,  0.10588f,  0.14902f, 0.79608f,   0.12549f,  0.15294f,
    0.81176f,  0.1451f,   0.15294f, 0.82745f,   0.16863f,  0.15294f,
    0.84314f,  0.18824f,  0.15294f, 0.85882f,   0.20784f,  0.15686f,
    0.87059f,  0.23137f,  0.16078f, 0.88627f,   0.2549f,   0.17255f,
    0.89804f,  0.27843f,  0.18039f, 0.9098f,    0.30196f,  0.19608f,
    0.92157f,  0.32941f,  0.20784f, 0.93333f,   0.35294f,  0.22353f,
    0.94118f,  0.37647f,  0.23529f, 0.94902f,   0.40392f,  0.25098f,
    0.95686f,  0.42745f,  0.26275f, 0.96078f,   0.45098f,  0.27451f,
    0.96863f,  0.47843f,  0.28627f, 0.97255f,   0.50588f,  0.29804f,
    0.97647f,  0.53333f,  0.30588f, 0.98039f,   0.55686f,  0.31765f,
    0.98431f,  0.58431f,  0.32941f, 0.98824f,   0.61176f,  0.34118f,
    0.98824f,  0.63529f,  0.35294f, 0.99216f,   0.65882f,  0.36863f,
    0.99216f,  0.68235f,  0.38039f, 0.99216f,   0.70588f,  0.39608f,
    0.99216f,  0.72549f,  0.40784f, 0.99216f,   0.74902f,  0.42353f,
    0.99608f,  0.76863f,  0.43922f, 0.99608f,   0.78824f,  0.45882f,
    0.99608f,  0.80784f,  0.47451f, 0.99608f,   0.82745f,  0.4902f,
    0.99608f,  0.84706f,  0.5098f,  0.99608f,   0.86275f,  0.52549f,
    0.99608f,  0.87843f,  0.5451f,  0.99608f,   0.89412f,  0.56471f,
    0.99608f,  0.9098f,   0.58824f, 0.99608f,   0.92549f,  0.61569f,
    0.99608f,  0.94118f,  0.64314f, 1.f,        0.95686f,  0.67059f,
    1.f,       0.97255f,  0.69412f, 1.f,        0.98431f,  0.71765f,
    1.f,       0.99216f,  0.73333f, 1.f,        0.99608f,  0.7451f,
    1.f,       1.f,       0.74902f, 0.99608f,   1.f,       0.7451f,
    0.98824f,  0.99608f,  0.72941f, 0.97255f,   0.98824f,  0.70588f,
    0.95686f,  0.98039f,  0.67843f, 0.93725f,   0.97255f,  0.65098f,
    0.91373f,  0.96471f,  0.61961f, 0.8902f,    0.95294f,  0.59216f,
    0.87059f,  0.9451f,   0.56471f, 0.85098f,   0.93725f,  0.5451f,
    0.83137f,  0.92941f,  0.52941f, 0.81569f,   0.92157f,  0.51373f,
    0.79608f,  0.91373f,  0.49412f, 0.77647f,   0.90588f,  0.47843f,
    0.75686f,  0.89804f,  0.46667f, 0.73725f,   0.88627f,  0.45098f,
    0.71373f,  0.87843f,  0.43922f, 0.69412f,   0.87059f,  0.43137f,
    0.67451f,  0.85882f,  0.41961f, 0.65098f,   0.85098f,  0.41569f,
    0.62745f,  0.84314f,  0.41176f, 0.60392f,   0.83137f,  0.40784f,
    0.58039f,  0.81961f,  0.40392f, 0.55686f,   0.81176f,  0.40392f,
    0.53333f,  0.8f,      0.4f,     0.50588f,   0.78824f,  0.4f,
    0.47843f,  0.77647f,  0.39608f, 0.4549f,    0.76471f,  0.39608f,
    0.42745f,  0.75294f,  0.39216f, 0.4f,       0.74118f,  0.38824f,
    0.37255f,  0.72941f,  0.38431f, 0.33725f,   0.71373f,  0.37647f,
    0.30588f,  0.70196f,  0.37255f, 0.27059f,   0.68627f,  0.36471f,
    0.23529f,  0.67451f,  0.35686f, 0.20392f,   0.65882f,  0.34902f,
    0.17255f,  0.64314f,  0.34118f, 0.1451f,    0.62745f,  0.32941f,
    0.12157f,  0.61176f,  0.32157f, 0.10196f,   0.59608f,  0.31373f,
    0.086275f, 0.58039f,  0.30588f, 0.070588f,  0.56078f,  0.29804f,
    0.058824f, 0.5451f,   0.28627f, 0.043137f,  0.52549f,  0.27843f,
    0.031373f, 0.50588f,  0.26667f, 0.023529f,  0.4902f,   0.25882f,
    0.011765f, 0.47059f,  0.24706f, 0.0078431f, 0.45098f,  0.23529f,
    0.f,       0.42745f,  0.22745f, 0.f,        0.40784f,  0.21569f};

ColorPalette ColorPalettes::MatlabHeatMap() {
  return std::move(
      ColorPalette(matlab_heat, sizeof(matlab_heat) / (sizeof(float) * 3)));
}

ColorPalette ColorPalettes::Batlow() {
  return std::move(ColorPalette(batlow, 99));
}

TrigonometricColorPalette ColorPalettes::Inigo6() {
  return std::move(TrigonometricColorPalette(
      {0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {2.0, 1.0, 0.0}, {0.50, 0.20, 0.25}));
}

ColorPalette::ColorPalette(const u32 *c, size_t n) {
  for (size_t i = 0; i < n; ++i)
    colors.emplace_back(Color::rbgFromU32(c[i]));
}

ColorPalette::ColorPalette(const f32 *c, size_t rgb_count) {
  for (size_t i = 0; i < rgb_count; i++) {
    float r = c[i * 3 + 0];
    float g = c[i * 3 + 1];
    float b = c[i * 3 + 2];
    colors.emplace_back(r, g, b);
  }
}

ColorPalette::ColorPalette() = default;

ColorPalette::ColorPalette(std::initializer_list<int> c) {
  size_t n = c.size() / 3;
  auto it = c.begin();
  for (size_t i = 0; i < n; i++) {
    float r = static_cast<float>(*it++) / 255.f;
    float g = static_cast<float>(*it++) / 255.f;
    float b = static_cast<float>(*it++) / 255.f;
    colors.emplace_back(r, g, b);
  }
}

ColorPalette::ColorPalette(std::initializer_list<double> c) : ColorPalette() {
  size_t n = c.size() / 3;
  auto it = c.begin();
  for (size_t i = 0; i < n; i++) {
    auto r = static_cast<float>(*it++);
    auto g = static_cast<float>(*it++);
    auto b = static_cast<float>(*it++);
    colors.emplace_back(r, g, b);
  }
}
} // namespace psa_anim
