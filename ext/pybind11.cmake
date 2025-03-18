include(FetchContent)

set(PYBIND11_PYTHON_VERSION 3.13.2 CACHE STRING "")
#find_package(
#  Python 3.8
#  COMPONENTS Interpreter Development
#  REQUIRED)

FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG        v2.13.6
)

FetchContent_MakeAvailable(pybind11)

