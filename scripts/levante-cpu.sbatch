#!/bin/bash

#SBATCH --account=ka1273
#SBATCH --job-name=scc_check
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --output=%x.%j.log
#SBATCH --time=00:30:00

ulimit -s unlimited
ulimit -c 0

#salloc --account=ka1273 --job-name=scc --partition=compute --ntasks=1 --time=00:10:00
# compiler flags
FLAGS='-O0 -nofma' 

# input file
IN_FILE=/work/ka1273/atm_R2B08.nc
# output file
OUT_FILE=/scratch/b/b383366/results/output.nc
# reference file associated with input file
REF_FILE=/scratch/b/b383366/results/atm_R2B08_ans_double.nc

# setup env
spack load cdo@2.5.1%gcc@=11.2.0
. $(pwd)/scripts/levante-setup.sh nvidia cpu
# build the code
BUILD=build_std_cpu_double
rm -rf $BUILD
cmake -B $BUILD -S . -DMU_IMPL=std -DMU_ARCH=x86_64 -DMU_ENABLE_MPI=ON -DMU_ENABLE_SINGLE=OFF -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_CXX_FLAGS="${FLAGS}"  && cmake --build $BUILD --parallel

# check correctness
. scripts/check-correctness.sh

# run the executable
BIN=$(pwd)/$BUILD/bin/graupel

srun -n 1 $BIN $IN_FILE $OUT_FILE
check $OUT_FILE $REF_FILE
#srun -n 1 $BIN $(pwd)/tasks/11k.nc output.nc
#check output.nc $(pwd)/reference_results/seq_11k_double.nc
#srun -n 1 $BIN $(pwd)/tasks/20k.nc output.nc
#check output.nc $(pwd)/reference_results/seq_20k_double.nc


#rm $OUT_FILE

echo "Job done!"
