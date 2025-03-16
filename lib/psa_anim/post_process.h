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
///\file post_process.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-12-07
///
///\brief

#ifndef psa_anim_psa_anim_SRC_POST_PROCESS_H
#define psa_anim_psa_anim_SRC_POST_PROCESS_H

#include <openfoam_poly_mesh.h>

namespace psa_anim {

/// \brief Smooths a scalar field based on limiting values
/// \param mesh
/// \param field_values
/// \param isovalue
std::vector<double> smooth(const PolyMesh *mesh,
                           const std::vector<double> &field_values,
                           double isovalue, size_t size);

} // namespace psa_anim

#endif // psa_anim_psa_anim_SRC_POST_PROCESS_H
