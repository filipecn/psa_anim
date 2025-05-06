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
///\file sim_path.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-03-22
///
///\brief

#ifndef PSA_ANIM_FOAM_UTILS_SIM_PATH_H
#define PSA_ANIM_FOAM_UTILS_SIM_PATH_H

// *********************************************************************************************************************
//                                                                                                            SimPath
// *********************************************************************************************************************
#include "hermes/logging/logging.h"
struct SimPath {
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  ///
  SimPath() = default;
  ///
  /// \param sim_folder
  explicit SimPath(const hermes::Path &sim_folder) { setRootPath(sim_folder); }
  // *******************************************************************************************************************
  //                                                                                                        OPERATORS
  // *******************************************************************************************************************
  ///
  void operator++() { current_frame_++; }
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  ///
  /// \return
  [[nodiscard]] bool finished() const {
    return current_frame_ >= frames_.size();
  }
  ///
  /// \return
  [[nodiscard]] const hermes::Path &currentFramePath() const {
    return frames_[current_frame_];
  }
  ///
  /// \return
  [[nodiscard]] size_t currentFrame() const { return current_frame_; }
  ///
  /// \return
  [[nodiscard]] const hermes::Path &root() const { return path_; }
  ///
  /// \return
  [[nodiscard]] const std::vector<hermes::Path> &frames() const {
    return frames_;
  }
  ///
  /// \param sim_folder
  void setRootPath(const hermes::Path &sim_folder) {
    path_ = sim_folder / "frames";
    auto folders = hermes::FileSystem::ls(
        path_, hermes::ls_options::sort | hermes::ls_options::directories);
    for (const auto &p : folders)
      if (hermes::Str::isNumber(p.name()) ||
          (hermes::Str::isInteger(p.name()) && std::stol(p.name()) != 0))
        frames_.emplace_back(p);
    if (!frames_.empty())
      std::sort(frames_.begin(), frames_.end(),
                [](const auto &a, const auto &b) -> bool {
                  return std::stof(a.name()) < std::stof(b.name());
                });
  }
  ///
  /// \param frame_id
  void setCurrentFrame(size_t frame_id) { current_frame_ = frame_id; }
  ///
  /// \param path
  void listFields(std::vector<std::string> &scalar_fields,
                  std::vector<std::string> &vector_fields,
                  const std::vector<std::string> &filter_names = {}) const {
    HERMES_LOG_AND_RETURN_IF_NOT(!frames_.empty(),
                                 "No simulation frames found!");
    auto cur_path = currentFramePath();
    if (!path_.isDirectory())
      return;
    auto files = hermes::FileSystem::ls(cur_path, hermes::ls_options::files);
    for (const auto &file : files) {
      // filter only
      if (file.hasExtension())
        continue;
      bool skip = true;
      for (const auto &name : filter_names)
        if (name == file.name())
          skip = false;
      //      if (file.name() != "fa_h" && file.name() != "fa_Us")
      //        continue;
      OpenFoamDict dict(file);
      auto field_type = dict["FoamFile"]["class"].value;
      if (field_type == "volVectorField")
        vector_fields.emplace_back(file.name());
      else if (field_type == "volScalarField")
        scalar_fields.emplace_back(file.name());
    }
    // HERMES_LOG("Listing {} fields", vector_fields.size() +
    // scalar_fields.size()); HERMES_LOG_ARRAY(scalar_fields);
    // HERMES_LOG_ARRAY(vector_fields);
  }

private:
  std::vector<hermes::Path> frames_{};
  hermes::Path path_{};
  size_t current_frame_{1};
};

#endif // PSA_ANIM_FOAM_UTILS_SIM_PATH_H
