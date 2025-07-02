#!/bin/bash

# compiler flags

# build
BUILD='build_std_cpu_single'

spack load /42ju4ng # netcdf-cxx4 compiled with nvhpc
module load nvhpc/24.7-gcc-11.2.0
module load openmpi/4.1.5-nvhpc-24.7
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64/


rm -rf $BUILD
cmake -B $BUILD -S . -DMU_IMPL=std -DMU_ARCH=x86_64 -DMU_ENABLE_SINGLE=ON -DMU_ENABLE_MPI=OFF -DCMAKE_CXX_COMPILER=nvc++  && cmake --build $BUILD --parallel