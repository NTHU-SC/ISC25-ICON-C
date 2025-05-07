#!/bin/bash

#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --output=test_v2.log
#SBATCH --time=01:00:00

ulimit -s unlimited
ulimit -c 0

#salloc --account=ka1273 --job-name=scc --partition=compute --ntasks=1 --time=00:10:00
# compiler flags
FLAGS='-O0 -stdpar=multicore -nofma' 

# input file
IN_FILE=/work/ka1273/atm_R2B08.nc
# output file
OUT_FILE=/scratch/b/b383366/results/output2.nc
# reference file associated with input file
REF_FILE=/scratch/b/b383366/results/atm_R2B08_ans_double.nc

# setup env
spack load cdo@2.5.1%gcc@=11.2.0
. $(pwd)/scripts/levante-setup.sh nvidia cpu

# build the code
PERCISION=single
BUILD=build_std_cpu_$PERCISION
. $(pwd)/scripts/build_cpu_$PERCISION.sh

# check correctness
. scripts/check-correctness.sh
export OMPI_MCA_opal_cuda_support=false

# run the executable
BIN=$(pwd)/$BUILD/bin/graupel

#srun -n 1 $BIN $IN_FILE $OUT_FILE
#check $OUT_FILE $REF_FILE
#for i in {1..16}; do
 #   echo "=== Run with -n $i ==="
  #  srun -n $i $BIN $(pwd)/tasks/11k.nc output.nc
   # check output.nc $(pwd)/reference_results/seq_11k_$PERCISION.nc
#
 #   srun -n $i $BIN $(pwd)/tasks/20k.nc output.nc
  #  check output.nc $(pwd)/reference_results/seq_20k_$PERCISION.nc
#
 #   srun -n $i $BIN $(pwd)/tasks/dbg.nc output.nc
  #  check output.nc $(pwd)/reference_results/seq_dbg_$PERCISION.nc
#done


#rm $OUT_FILE

echo "Job done!"
