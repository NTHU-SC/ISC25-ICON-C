#!/bin/bash

# load the default C++ compiler
module load gcc

# load netcdf support
spack load netcdf-cxx4@4.3.1

# load the nvidia libs & export path for C++ compiler
module load nvhpc/24.7-gcc-11.2.0
export LD_LIBRARY_PATH=/sw/spack-levante/llvm-18.1.6-br53hv/lib:/sw/spack-levante/gcc-13.3.0-s2dxrt/lib64/:$LD_LIBRARY_PATH

# load tool for comparing results
module load cdo