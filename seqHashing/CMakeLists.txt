cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(vtk-bs-hashing)

set(CMAKE_CXX_STANDARD 14)

set(FETCHCONTENT_QUIET OFF)

include(FetchContent)

FetchContent_Declare(
  VTK
  GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git
  GIT_TAG 15527a6e8ab744da0a52622fe7ab34ba6f823d8e
  GIT_PROGRESS ON
  SOURCE_DIR VTK
)
FetchContent_MakeAvailable(VTK)

add_executable(seqHashing main.cxx)
target_link_libraries(seqHashing PUBLIC VTK::CommonCore VTK::IOXML VTK::FiltersGeometry)

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS seqHashing
  MODULES VTK::CommonCore VTK::IOXML VTK::FiltersGeometry
)
