#!/bin/bash
 
c1="gnu"
c2="intel"
c3="nvidia"

a1="cpu"
a2="gpu"

module purge 
#spack load netcdf-cxx4@4.3.1

echo "compilation for architecture [$2] using compilation chain [$1]"
echo ""

if [[ "$1" == "$c1" ]]; then
    module load gcc
elif [[ "$1" == "$c2" ]]; then
    module load intel-oneapi-compilers/2022.0.1-gcc-11.2.0
elif [[ "$1" == "$c3" ]]; then
    spack load /42ju4ng # netcdf-cxx4 compiled with nvhpc
    module load nvhpc/24.7-gcc-11.2.0
    module load openmpi/4.1.5-nvhpc-24.7
    export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64/
    #export LD_LIBRARY_PATH=/sw/spack-levante/llvm-18.1.6-br53hv/lib:/sw/spack-levante/gcc-13.3.0-s2dxrt/lib64/:$LD_LIBRARY_PATH 
else 
    echo "$1 is not supported!"
fi

if [[ "$2" == "$a2" ]]; then 
    module load nvhpc/24.7-gcc-11.2.0
    export LD_LIBRARY_PATH=/sw/spack-levante/llvm-18.1.6-br53hv/lib:/sw/spack-levante/gcc-13.3.0-s2dxrt/lib64/:$LD_LIBRARY_PATH 
fi
