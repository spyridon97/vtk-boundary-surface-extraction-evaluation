#!/bin/bash

#set -x
set -e

base_dir=`pwd`

module load git
module load ninja
module load cmake
module load gcc-native
module load rocm

kokkos_repo=https://github.com/kokkos/kokkos.git
#kokkos_branch=3.7.02
kokkos_branch=4.4.01

#cmake_generator="Unix Makefiles"
cmake_generator="Ninja"

kokkos_src_dir=${base_dir}/kokkos/src
kokkos_build_dir=${base_dir}/kokkos/build
kokkos_install_dir=${base_dir}/kokkos/install

if true; then
rm -rf "${kokkos_src_dir}"
rm -rf "${kokkos_build_dir}"
rm -rf "${kokkos_install_dir}"
if [[ ! -d ${kokkos_src_dir} ]] ; then
  git clone -b ${kokkos_branch} ${kokkos_repo} "${kokkos_src_dir}"
fi
cmake -G "${cmake_generator}" -S "${kokkos_src_dir}" -B "${kokkos_build_dir}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON\
  -DKokkos_ARCH_VEGA90A=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=OFF \
  -DCMAKE_INSTALL_PREFIX="${kokkos_install_dir}" \
  -DCMAKE_CXX_FLAGS="--amdgpu-target=gfx90a" \
  -DCMAKE_CXX_COMPILER=hipcc
time cmake --build "${kokkos_build_dir}" -j
cmake --install "${kokkos_build_dir}"
fi

tbb_repo=https://github.com/oneapi-src/oneTBB.git
tbb_branch=v2021.13.0

tbb_src_dir=${base_dir}/tbb/src
tbb_build_dir=${base_dir}/tbb/build
tbb_install_dir=${base_dir}/tbb/install

if true; then
rm -rf "${tbb_src_dir}"
rm -rf "${tbb_build_dir}"
rm -rf "${tbb_install_dir}"
if [[ ! -d ${tbb_src_dir} ]] ; then
  git clone -b ${tbb_branch} ${tbb_repo} "${tbb_src_dir}"
fi
cmake -G "${cmake_generator}" -S "${tbb_src_dir}" -B "${tbb_build_dir}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON\
  -DTBB_TEST=OFF\
  -DCMAKE_INSTALL_PREFIX="${tbb_install_dir}" \
  -DCMAKE_CXX_COMPILER=amdclang++ \
  -DCMAKE_C_COMPILER=amdclang
time cmake --build "${tbb_build_dir}" -j
cmake --install "${tbb_build_dir}"
fi

vtkefe_src_dir=${base_dir}
vtkefe_build_dir=${base_dir}/build
vtkefe_install_dir=${base_dir}/install

if true; then
rm -rf "${vtkefe_build_dir}"
rm -rf "${vtkefe_install_dir}"
cmake -G "${cmake_generator}" -S "${vtkefe_src_dir}" -B "${vtkefe_build_dir}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON\
  -DVTKm_ENABLE_KOKKOS=ON \
  -DVTKm_USE_DOUBLE_PRECISION=ON \
  -DCMAKE_HIP_ARCHITECTURES="gfx90a" \
  -DCMAKE_PREFIX_PATH="${kokkos_install_dir}" \
  -DTBB_ROOT="${tbb_install_dir}" \
  -DTBB_INCLUDE_DIR="${tbb_install_dir}/include" \
  -DTBB_LIBRARY="${tbb_install_dir}/lib64/libtbb.so" \
  -DCMAKE_INSTALL_PREFIX="${vtkefe_install_dir}" \
  -DCMAKE_CXX_COMPILER=amdclang++ \
  -DCMAKE_C_COMPILER=amdclang

time cmake --build "${vtkefe_build_dir}" -j
cmake --install "${vtkefe_build_dir}"
fi
