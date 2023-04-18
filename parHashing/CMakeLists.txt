cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(vtk-par-bs-hashing)

set(CMAKE_CXX_STANDARD 14)

set(FETCHCONTENT_QUIET OFF)

include(FetchContent)

FetchContent_Declare(
  VTK
  GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git
  GIT_TAG b2bdf26f6b45b38b08f105374306664bd0a1843c
  GIT_PROGRESS ON
  SOURCE_DIR VTK
)
FetchContent_MakeAvailable(VTK)

add_executable(parHashing main.cxx)
target_link_libraries(parHashing PUBLIC VTK::CommonCore VTK::IOXML VTK::FiltersGeometry)

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS parHashing
  MODULES VTK::CommonCore VTK::IOXML VTK::FiltersGeometry
)
