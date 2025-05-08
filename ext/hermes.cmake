include(FetchContent)

FetchContent_Declare(
  hermes
  GIT_REPOSITORY https://github.com/filipecn/hermes.git
  GIT_TAG "v1.0.0"
)

FetchContent_MakeAvailable(hermes)

SET(HERMES_INCLUDES ${hermes_SOURCE_DIR})
