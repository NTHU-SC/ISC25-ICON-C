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
FLAGS='-O0 -nofma' 

# input file
IN_FILE=$(pwd)/tasks/dbg.nc
# output file
OUT_FILE=output_dbg_double.nc
# reference file associated with input file
REF_FILE=$(pwd)/reference_results/seq_dbg_double.nc

# setup env
. scripts/levante-setup.sh 

# build the code
cmake -B build_std_cpu -S . -DMU_IMPL=std -DMU_ARCH=x86_64 -DMU_ENABLE_SINGLE=OFF -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_CXX_FLAGS="${FLAGS}"  && cmake --build build_std_cpu --parallel

# cleanup previous results (if any)
if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
fi 

# run the executable
./build_std_cpu/bin/graupel $IN_FILE $OUT_FILE

# check correctness
. scripts/check-correctness.sh
check $OUT_FILE $REF_FILE

echo "Job done!"
