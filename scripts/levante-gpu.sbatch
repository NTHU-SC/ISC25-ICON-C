#!/bin/bash
  
#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --output=%x.%j.log
#SBATCH --exclusive
#SBATCH --time=00:10:00

# compiler flags
FLAGS='-O0 -stdpar=gpu -gpu=cc80,nofma,noflushz' 

# input file
IN_FILE=$(pwd)/tasks/dbg.nc
# output file
OUT_FILE=output_gpu_double.nc
# reference file associated with input file
REF_FILE=$(pwd)/reference_results/seq_dbg_double.nc

# setup env
. scripts/levante-setup.sh 

# build the code
cmake -B build_std_gpu -S . -DMU_IMPL=std -DMU_ARCH=a100 -DMU_ENABLE_SINGLE=OFF -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_CXX_FLAGS="${FLAGS}" && cmake --build build_std_gpu --parallel

# cleanup previous results (if any)
if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
fi 

# run the executable
./build_std_gpu/bin/graupel $IN_FILE $OUT_FILE

# check correctness
. scripts/check-correctness.sh
check $OUT_FILE $REF_FILE

echo "Job done!"
