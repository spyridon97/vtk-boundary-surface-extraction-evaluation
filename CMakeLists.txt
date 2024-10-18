cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(vtk-external-facelist-evaluation)

# VTK flags
set(VTK_SMP_IMPLEMENTATION_TYPE TBB CACHE STRING "" FORCE)
#set(VTK_USE_KOKKOS ON CACHE BOOL "" FORCE)
#set(VTK_KOKKOS_BACKEND HIP CACHE STRING "" FORCE)
set(VTK_MODULE_ENABLE_VTK_vtkvtkm YES CACHE BOOL "" FORCE)
set(VTK_MODULE_ENABLE_VTK_AcceleratorsVTKmCore YES CACHE BOOL "" FORCE)
set(VTK_MODULE_ENABLE_VTK_AcceleratorsVTKmDataModel YES CACHE BOOL "" FORCE)
set(VTK_MODULE_ENABLE_VTK_AcceleratorsVTKmFilters YES CACHE BOOL "" FORCE)

set(FETCHCONTENT_QUIET OFF)
include(FetchContent)

FetchContent_Declare(
  VTK
  GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git
  GIT_TAG v9.4.0.rc2 # Includes vtkm 2.2 and it's VTK 9.3.1+
  GIT_PROGRESS ON
  SOURCE_DIR VTK
)
FetchContent_MakeAvailable(VTK)

FetchContent_Declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.4.2
  GIT_PROGRESS ON
  SOURCE_DIR CLI11
)
FetchContent_MakeAvailable(CLI11)

set(headers
  src/CellFaceMinPointId.h

  src/ExternalFacesHashCountFnv1a.h
  src/ExternalFacesHashCountMinPointId.h
  src/ExternalFacesHashFightFnv1a.h
  src/ExternalFacesHashFightMinPointId.h
  src/ExternalFacesHashSortFnv1a.h
  src/ExternalFacesHashSortMinPointId.h

  src/vtkGeometryFilterSClassifier.h
  src/vtkDataSetSurfaceFilterSHash.h
  src/vtkGeometryFilterPClassifier.h
  src/vtkDataSetSurfaceFilterSHash.h

  src/YamlWriter.h
  src/Arguments.h
)

set(sources
  src/vtkGeometryFilterSClassifier.cxx
  src/vtkDataSetSurfaceFilterSHash.cxx
  src/vtkGeometryFilterPClassifier.cxx
  src/vtkGeometryFilterPHash.cxx

  src/Arguments.cxx

  src/ExternalFacelistEvaluation.cxx
)

add_executable(${PROJECT_NAME} ${headers} ${sources})

target_include_directories(${PROJECT_NAME} PRIVATE ${VTKm_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC
  VTK::CommonCore VTK::IOXML VTK::FiltersGeometry VTK::vtkvtkm VTK::AcceleratorsVTKmDataModel CLI11::CLI11)

set_source_files_properties(${headers}
  PROPERTIES HEADER_FILE_ONLY TRUE
)
