#!/bin/bash

# compiler flags
FLAGS='-O3 -march=native -stdpar=gpu -gpu=cc80,nofma' 

# build
BUILD='build_std_gpu_mpi_double'

spack load /42ju4ng # netcdf-cxx4 compiled with nvhpc
module load nvhpc/24.7-gcc-11.2.0
module load openmpi/4.1.5-nvhpc-24.7
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64/

module load nvhpc/24.7-gcc-11.2.0
export LD_LIBRARY_PATH=/sw/spack-levante/llvm-18.1.6-br53hv/lib:/sw/spack-levante/gcc-13.3.0-s2dxrt/lib64/:$LD_LIBRARY_PATH 

rm -rf $BUILD

# build the code
cmake -B $BUILD -S . -DMU_IMPL=std -DMU_ENABLE_MPI=ON -DMU_ARCH=a100 -DMU_ENABLE_SINGLE=OFF -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_CXX_FLAGS="${FLAGS}" && cmake --build $BUILD --parallel