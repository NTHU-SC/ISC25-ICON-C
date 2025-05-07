#!/bin/bash

#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --output=%x.%j.log
#SBATCH --time=00:10:00

ulimit -s unlimited
ulimit -c 0

# compiler flags
FLAGS='-O0 -stdpar=multicore -nofma'

# build
BUILD='build_std_cpu_double'

spack load /42ju4ng # netcdf-cxx4 compiled with nvhpc
module load nvhpc/24.7-gcc-11.2.0
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64/


rm -rf $BUILD
cmake -B $BUILD -S . -DMU_IMPL=std -DMU_ARCH=x86_64 -DMU_ENABLE_SINGLE=OFF -DMU_ENABLE_MPI=ON -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_CXX_FLAGS="${FLAGS}"  && cmake --build $BUILD --parallel

