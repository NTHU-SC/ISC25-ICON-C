#!/bin/bash
  
#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --output=atm_v7.log
#SBATCH --exclusive
#SBATCH --time=00:30:00

# compiler flags
FLAGS='-O3 -march=native -stdpar=gpu -gpu=cc80,nofma' 

# input file
IN_FILE=/work/ka1273/atm_R2B08.nc
# output file
OUT_FILE=/scratch/b/b383366/outputs/output.nc
# reference file associated with input file
REF_FILE=/scratch/b/b383366/results/atm_R2B08_ans_double.nc

BUILD=build_std_gpu
EXE=$(pwd)/$BUILD/bin/graupel

# setup env
spack load cdo@2.5.1%gcc@=11.2.0
spack load /42ju4ng # netcdf-cxx4 compiled with nvhpc
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64/
module load nvhpc/24.7-gcc-11.2.0
export LD_LIBRARY_PATH=/sw/spack-levante/llvm-18.1.6-br53hv/lib:/sw/spack-levante/gcc-13.3.0-s2dxrt/lib64/:$LD_LIBRARY_PATH


# build the code
rm -rf $BUILD
cmake -B $BUILD -S . -DMU_IMPL=std -DMU_ARCH=a100 -DMU_ENABLE_MPI=ON -DMU_ENABLE_SINGLE=OFF -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_CXX_FLAGS="${FLAGS}" && cmake --build $BUILD --parallel

# cleanup previous results (if any)
if [ -f $OUT_FILE ]; then
    rm $OUT_FILE
fi 


#. scripts/check-correctness.sh
#srun -n 4 $EXE $(pwd)/tasks/11k.nc $OUT_FILE
#check $OUT_FILE $(pwd)/reference_results/seq_11k_double.nc
#srun -n 4 $EXE $(pwd)/tasks/20k.nc $OUT_FILE
#check $OUT_FILE $(pwd)/reference_results/seq_20k_double.nc
#srun -n 4 $EXE $(pwd)/tasks/dbg.nc $OUT_FILE
#check $OUT_FILE $(pwd)/reference_results/seq_dbg_double.nc
export MULTI_GRAUPEL="1"
IN_DIR=$BUILD/bin

: ${no_of_nodes:=${SLURM_JOB_NUM_NODES:=1}} ${mpi_procs_pernode:=1}
export no_of_nodes
num_io_procs=0
((mpi_total_procs=no_of_nodes * mpi_procs_pernode))

srun -l --kill-on-bad-exit=1 --nodes=${SLURM_JOB_NUM_NODES:-1} --gpus-per-node=1 --ntasks=$((no_of_nodes * mpi_procs_pernode)) $(pwd)/wrapper.sh $IN_DIR/graupel $IN_FILE $OUT_FILE

. scripts/check-correctness.sh
check $OUT_FILE $REF_FILE

echo "Job done!"
