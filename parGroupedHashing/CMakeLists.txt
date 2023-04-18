cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(vtk-par-bs-grouped-hashing)

set(CMAKE_CXX_STANDARD 14)

set(FETCHCONTENT_QUIET OFF)

include(FetchContent)

FetchContent_Declare(
  VTK
  GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git
  GIT_TAG 38f3080453c760acfd5e9c0c9dfff18e1b177e4e
  GIT_PROGRESS ON
  SOURCE_DIR VTK
)
FetchContent_MakeAvailable(VTK)

add_executable(parGroupedHashing main.cxx)
target_link_libraries(parGroupedHashing PUBLIC VTK::CommonCore VTK::IOXML VTK::FiltersGeometry)

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS parGroupedHashing
  MODULES VTK::CommonCore VTK::IOXML VTK::FiltersGeometry
)
