cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(vtk-par-bs-classifier)

set(CMAKE_CXX_STANDARD 14)

set(FETCHCONTENT_QUIET OFF)

include(FetchContent)

FetchContent_Declare(
  VTK
  GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git
  GIT_TAG c44fd2ba0d86779d65f5eda56a1aaceff2954807
  GIT_PROGRESS ON
  SOURCE_DIR VTK
)
FetchContent_MakeAvailable(VTK)

add_executable(parClassifier main.cxx)
target_link_libraries(parClassifier PUBLIC VTK::CommonCore VTK::IOXML VTK::FiltersGeometry)

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS parClassifier
  MODULES VTK::CommonCore VTK::IOXML VTK::FiltersGeometry
)
