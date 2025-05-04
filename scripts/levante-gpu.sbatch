#!/bin/bash
  
#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --output=%x.%j.log
#SBATCH --exclusive
#SBATCH --time=00:10:00

# compiler flags
FLAGS='-O3 -stdpar=gpu' 

# input file
IN_FILE=$(pwd)/tasks/dbg.nc
# output file
OUT_FILE=output_gpu_double.nc
# reference file associated with input file
REF_FILE=$(pwd)/reference_results/seq_dbg_double.nc

EXE=$(pwd)/build_std_gpu/bin/graupe

# setup env
. scripts/levante-setup.sh nvidia gpu

# build the code
rm -rf build_std_gpu
cmake -B build_std_gpu -S . -DMU_IMPL=std -DMU_ARCH=a100 -DMU_ENABLE_SINGLE=OFF -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_CXX_FLAGS="${FLAGS}" && cmake --build build_std_gpu --parallel

# cleanup previous results (if any)
if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
fi 


. scripts/check-correctness.sh
$EXE $IN_FILE $OUT_FILE
check $OUT_FILE $REF_FILE

echo "Job done!"
