#!/bin/bash

build_type=DEBUG
eamxx_src=$HOME/scream-mam4/components/eamxx

#
# copy to your eamxx build directory, e.g., ~/scream/components/eamxx/build
#
cmake -Wno-dev \
  -D CMAKE_BUILD_TYPE=$build_type \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_Fortran_COMPILER=mpifort \
  -D CMAKE_Fortran_FLAGS="-fallow-argument-mismatch" \
  -D CMAKE_CXX_STANDARD=17 \
  -D MPIEXEC_EXECUTABLE=`which mpiexec` \
  -D SCREAM_CIME_BUILD=OFF \
  -D SCREAM_DOUBLE_PRECISION:BOOL=TRUE \
  -D SCREAM_INPUT_ROOT:PATH=/data/scream-input \
  -D SCREAM_DYNAMICS_DYCORE=HOMME \
  -D SCREAM_ENABLE_MAM=ON \
  -D Kokkos_ENABLE_LIBDL=OFF \
  -D Kokkos_ENABLE_DEBUG=TRUE \
  -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=OFF \
  -D Kokkos_ENABLE_SERIAL=ON \
  -D Kokkos_ENABLE_OPENMP=ON \
  -D Kokkos_ENABLE_DEPRECATED_CODE=OFF \
  -D NetCDF_C_PATHS=$NETCDF_ROOT \
  -D NetCDF_Fortran_PATHS=$NETCDF_ROOT \
  -D PnetCDF_C_PATHS=$PNETCDF_ROOT \
  -D PnetCDF_Fortran_PATHS=$PNETCDF_ROOT \
$eamxx_src
