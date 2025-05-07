#!/bin/bash
  
#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --output=gpu_latest.log
#SBATCH --exclusive
#SBATCH --time=00:10:00

# compiler flags
FLAGS='-O3 -march=native -stdpar=gpu -gpu=cc80,nofma' 

# input file
IN_FILE=/work/ka1273/atm_R2B08.nc
# output file
OUT_FILE=/scratch/b/b383366/outputs/output.nc
# reference file associated with input file
REF_FILE=/scratch/b/b383366/results/atm_R2B08_ans_double.nc

EXE=$(pwd)/$BUILD/bin/graupel

# setup env
spack load cdo@2.5.1%gcc@=11.2.0
spack load /42ju4ng # netcdf-cxx4 compiled with nvhpc
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64/
module load nvhpc/24.7-gcc-11.2.0
export LD_LIBRARY_PATH=/sw/spack-levante/llvm-18.1.6-br53hv/lib:/sw/spack-levante/gcc-13.3.0-s2dxrt/lib64/:$LD_LIBRARY_PATH


# build the code
PERCISION=double
BUILD=build_std_gpu_$PERCISION
. $(pwd)/scripts/build_gpu_$PERCISION.sh

# cleanup previous results (if any)
if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
fi 

export MULTI_GRAUPEL="1"
IN_DIR=$BUILD/bin

. scripts/check-correctness.sh


srun -l --kill-on-bad-exit=1 --nodes=1 --gpus-per-node=1 --ntasks=1 $(pwd)/wrapper.sh $IN_DIR/graupel $IN_FILE $OUT_FILE
check $OUT_FILE $REF_FILE

#srun -l --kill-on-bad-exit=1 --nodes=1 --gpus-per-node=1 --ntasks=1 $IN_DIR/graupel $(pwd)/tasks/11k.nc $OUT_FILE
#check $OUT_FILE $(pwd)/reference_results/seq_11k_$PERCISION.nc


#srun -l --kill-on-bad-exit=1 --nodes=1 --gpus-per-node=1 --ntasks=1 $IN_DIR/graupel $(pwd)/tasks/20k.nc $OUT_FILE
#check $OUT_FILE $(pwd)/reference_results/seq_20k_$PERCISION.nc

#srun -l --kill-on-bad-exit=1 --nodes=1 --gpus-per-node=1 --ntasks=1 $IN_DIR/graupel $(pwd)/tasks/dbg.nc $OUT_FILE
#check $OUT_FILE $(pwd)/reference_results/seq_dbg_$PERCISION.nc

#. scripts/check-correctness.sh
#check $OUT_FILE $REF_FILE

echo "Job done!"
