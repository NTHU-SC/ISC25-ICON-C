#! /bin/bash

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#
# Based on the run_wrapper for JUWELS Booster by Luis Kornblueh
# Adapted to Levante (only 2 NICs instead of 4 as on JUWELS Booster)
#
# nvidia-smi topo -mp
#	GPU0	GPU1	GPU2	GPU3	mlx5_0	mlx5_1	CPU Affinity	NUMA Affinity
# GPU0	 X 	SYS	SYS	SYS	SYS	SYS	48-63,176-191	3
# GPU1	SYS	 X 	SYS	SYS	PIX	SYS	16-31,144-159	1
# GPU2	SYS	SYS	 X 	SYS	SYS	PIX	112-127,240-255	7
# GPU3	SYS	SYS	SYS	 X 	SYS	SYS	80-95,208-223	5
# mlx5_0	SYS	PIX	SYS	SYS	 X 	SYS		
# mlx5_1	SYS	SYS	PIX	SYS	SYS	 X 		
#
# Legend:
#
#  X    = Self
#  SYS  = Connection traversing PCIe as well as the SMP interconnect between NUMA nodes (e.g., QPI/UPI)
#  NODE = Connection traversing PCIe as well as the interconnect between PCIe Host Bridges within a NUMA node
#  PHB  = Connection traversing PCIe as well as a PCIe Host Bridge (typically the CPU)
#  PXB  = Connection traversing multiple PCIe bridges (without traversing the PCIe Host Bridge)
#  PIX  = Connection traversing at most a single PCIe bridge
#
#____________________________________________________________________________________________________


while getopts n:o:e:f:g: argv
do
    case "${argv}" in
        n) mpi_total_procs=${OPTARG};;
        o) io_tasks=${OPTARG};;
        e) executable=${OPTARG};;
        f) input_file=${OPTARG};;
        g) output_file=${OPTARG};;
    esac
done

set -eu

lrank=$SLURM_LOCALID%4

export OMPI_MCA_pml=ucx
export OMPI_MCA_btl="^vader,tcp,openib,smcuda"

# need to check in run script that the variables make sense and are
# exported!

(( compute_tasks = mpi_total_procs - io_tasks ))

if (( SLURM_PROCID < compute_tasks ))
then

    echo Compute process $SLURM_LOCALID on $(hostname)

    numanode=(2-3 0-1 6-7 4-5)
    gpus=(0 1 2 3)
    nics=(mlx5_0:1 mlx5_0:1 mlx5_1:1 mlx5_1:1)
    reorder=(0 1 2 3)

    nic_reorder=(${nics[${reorder[0]}]}
                 ${nics[${reorder[1]}]}
                 ${nics[${reorder[2]}]}
                 ${nics[${reorder[3]}]})
    numanode_reorder=(${numanode[${reorder[0]}]}
                      ${numanode[${reorder[1]}]}
                      ${numanode[${reorder[2]}]}
                      ${numanode[${reorder[3]}]})

    export UCX_NET_DEVICES=${nic_reorder[lrank]}
    export CUDA_VISIBLE_DEVICES=${gpus[${reorder[lrank]}]}

    export UCX_RNDV_SCHEME=put_zcopy
    export UCX_RNDV_THRESH=16384

    export UCX_IB_GPU_DIRECT_RDMA=yes

    export UCX_TLS=all
    export UCX_MEMTYPE_CACHE=n

else

    echo IO process $SLURM_LOCALID on $(hostname)

    numanode=(2-3 0-1 6-7 4-5)
    nics=(mlx5_0:1 mlx5_0:1 mlx5_1:1 mlx5_1:1)    
    reorder=(0 1 2 3)
 
    nic_reorder=(${nics[${reorder[0]}]}
                 ${nics[${reorder[1]}]}
                 ${nics[${reorder[2]}]}
                 ${nics[${reorder[3]}]}) 

    numanode_reorder=(${numanode[${reorder[0]}]}
                      ${numanode[${reorder[1]}]}
                      ${numanode[${reorder[2]}]}
                      ${numanode[${reorder[3]}]}) 
    
    export UCX_NET_DEVICES=${nic_reorder[lrank]}

    export UCX_RNDV_SCHEME=put_zcopy
    export UCX_RNDV_THRESH=16384

    export UCX_IB_GPU_DIRECT_RDMA=yes

    export UCX_TLS=cma,rc,mm,cuda_ipc,cuda_copy,gdr_copy    
    export UCX_MEMTYPE_CACHE=n
    
fi


#----------------------------------------------------------- nvsmi-logger --------------------------------------------
lrankb=$SLURM_LOCALID
nvsmi_logger_PID=0
function kill_nvsmi()
{
    set +x
    if (( nvsmi_logger_PID != 0 ))
    then
        kill $nvsmi_logger_PID
    fi
}
trap kill_nvsmi ERR

# Local Task 0 runs always the nvidia-smi logger
if (( lrankb == 0 ))
then
    set +x
    logdir=${logdir=.}
    # Start logger in background. It will be killed by the ERR trap or kill_nvsmi.
    loop_repetition_time=500000000 # in nano seconds
    #while sleep 0.$(( ( 1999999999 - 1$(date +%N) ) % loop_repetition_time ))
    while true
    do
        LC_TIME=en_US date -Ins
        nvidia-smi --format=csv --query-gpu=index,power.draw,utilization.gpu,temperature.gpu,memory.used
    done > $logdir/nvsmi.log.$SLURM_JOB_ID.$SLURM_NODEID &
    nvsmi_logger_PID=$!
    set -x
fi

#----------------------------------------------------------- nvsmi-logger --------------------------------------------

numactl --cpunodebind=${numanode_reorder[$lrank]} --membind=${numanode_reorder[$lrank]} $executable $input_file $output_file

kill_nvsmi
