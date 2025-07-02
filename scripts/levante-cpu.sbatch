#!/bin/bash

#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --output=cpu_double_check.log
#SBATCH --time=00:30:00

ulimit -s unlimited
ulimit -c 0

export OMPI_MCA_opal_cuda_support=false


PERCISION=double

# setup env
#. scripts/levante-setup.sh nvidia cpu

# build the code

#BUILD=build_std

IN_FILE=/work/ka1273/atm_R2B08.nc

OUT_FILE=/scratch/b/b383366/outputs/output_cpu_$PERCISION.nc

REF_FILE=/scratch/b/b383366/results/atm_R2B08_ans_$PERCISION.nc


# build
BUILD=build_std_cpu_$PERCISION

spack load /42ju4ng # netcdf-cxx4 compiled with nvhpc
module load nvhpc/24.7-gcc-11.2.0
module load openmpi/4.1.5-nvhpc-24.7
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64/

# build the code
. scripts/build_cpu_$PERCISION.sh

module load cdo

# cleanup previous results (if any)
if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
fi 

. scripts/check-correctness.sh

# run the executable
# BIN=$(pwd)/build/bin/graupel

BIN=$(pwd)/$BUILD/bin/graupel

# $BIN $IN_FILE $OUT_FILE
# check $OUT_FILE $REF_FILE

for NP in $(seq 1 1); do
  echo "=== Running with $NP MPI task(s) ==="

  srun --ntasks=$NP $BIN $(pwd)/tasks/dbg.nc $OUT_FILE
  #$BIN $(pwd)/tasks/dbg.nc $(pwd)/seq/seq_dbg_$PERCISION.nc
  check $OUT_FILE /home/b/b383366/clean-repo/seq/seq_dbg_$PERCISION.nc

  rm -rf $OUT_FILE

  srun --ntasks=$NP $BIN $(pwd)/tasks/11k.nc $OUT_FILE
  #$BIN $(pwd)/tasks/11k.nc $(pwd)/seq/seq_11k_$PERCISION.nc
  check $OUT_FILE /home/b/b383366/clean-repo/seq/seq_11k_$PERCISION.nc

  rm -rf $OUT_FILE

  srun --ntasks=$NP $BIN $(pwd)/tasks/20k.nc $OUT_FILE
  #$BIN $(pwd)/tasks/20k.nc $(pwd)/seq/seq_20k_$PERCISION.nc
  check $OUT_FILE /home/b/b383366/clean-repo/seq/seq_20k_$PERCISION.nc

  rm -rf $OUT_FILE

  echo ""
done

srun --ntasks=1 $BIN $IN_FILE $OUT_FILE
check $OUT_FILE $REF_FILE

echo "All runs completed."
