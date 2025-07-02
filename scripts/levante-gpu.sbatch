#!/bin/bash
  
#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --output=%x.%j.log
#SBATCH --exclusive
#SBATCH --time=00:10:00

PERCISION=single

IN_FILE=/work/ka1273/atm_R2B08.nc

OUT_FILE=/scratch/b/b383366/outputs/output_gpu_$PERCISION.nc

REF_FILE=/scratch/b/b383366/results/atm_R2B08_ans_$PERCISION.nc

# setup env
. scripts/levante-setup.sh nvidia gpu

# build the code
. scripts/build_gpu_$PERCISION.sh

module load cdo

# cleanup previous results (if any)
if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
fi 

. scripts/check-correctness.sh

# run the executable

BIN=$(pwd)/$BUILD/bin/graupel

srun -l --kill-on-bad-exit=1 $BIN $(pwd)/tasks/dbg.nc $OUT_FILE
check $OUT_FILE $(pwd)/reference_results/seq_dbg_$PERCISION.nc

srun -l --kill-on-bad-exit=1 $BIN $(pwd)/tasks/11k.nc $OUT_FILE
check $OUT_FILE  $(pwd)/reference_results/seq_11k_$PERCISION.nc

srun -l --kill-on-bad-exit=1 $BIN $(pwd)/tasks/20k.nc $OUT_FILE
check $OUT_FILE  $(pwd)/reference_results/seq_20k_$PERCISION.nc

srun -l --kill-on-bad-exit=1 $BIN $IN_FILE $OUT_FILE
check $OUT_FILE $REF_FILE

echo "Job done!"
