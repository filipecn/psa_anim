include(FetchContent)

FetchContent_Declare(
  hermes
  GIT_REPOSITORY https://github.com/filipecn/hermes.git
  GIT_TAG main
  #GIT_TAG "v.1.0.0-beta"
)

FetchContent_MakeAvailable(hermes)

SET(HERMES_INCLUDES ${hermes_SOURCE_DIR})
