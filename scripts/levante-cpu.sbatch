#!/bin/bash

#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=compute
#SBATCH --ntasks=3
#SBATCH --output=%x.%j.log
#SBATCH --time=00:10:00

ulimit -s unlimited
ulimit -c 0

# compiler flags
FLAGS='-O0 -nofma' 

# input file
IN_FILE=$(pwd)/tasks/dbg.nc
# output file
OUT_FILE=output_dbg_double.nc
# reference file associated with input file
REF_FILE=$(pwd)/reference_results/seq_dbg_single.nc

# setup env
ml purge
spack load /42ju4ng # netcdf-cxx4 compiled with nvhpc
module load nvhpc/24.7-gcc-11.2.0
module load openmpi/4.1.5-nvhpc-24.7
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64/
module load nvhpc/24.7-gcc-11.2.0
export LD_LIBRARY_PATH=/sw/spack-levante/llvm-18.1.6-br53hv/lib:/sw/spack-levante/gcc-13.3.0-s2dxrt/lib64/:$LD_LIBRARY_PATH
spack load cdo@2.5.1%gcc@=11.2.0

# build the code
rm -rf build_std_cpu
cmake -B build_std_cpu -S . -DMU_IMPL=std -DMU_ARCH=x86_64 -DMU_ENABLE_MPI=ON -DMU_ENABLE_SINGLE=ON -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_CXX_FLAGS="${FLAGS}"  && cmake --build build_std_cpu --parallel

# cleanup previous results (if any)
if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
fi 

# check correctness
. scripts/check-correctness.sh

# run the executable
BIN=$(pwd)/build_std_cpu/bin/graupel
$BIN $IN_FILE $OUT_FILE
check $OUT_FILE $REF_FILE
$BIN $(pwd)/tasks/11k.nc output.nc
check output.nc $(pwd)/reference_results/seq_11k_single.nc
$BIN $(pwd)/tasks/20k.nc output.nc
check output.nc $(pwd)/reference_results/seq_20k_single.nc


rm $OUT_FILE

echo "Job done!"
