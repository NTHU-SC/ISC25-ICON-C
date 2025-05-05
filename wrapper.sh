#! /bin/bash

# Output to ./nsys_reports/rank_$N.nsys-rep
nsys  profile \
	-o "./nsys_reports/gpu_$SLURM_PROCID.nsys-rep" \
	--mpi-impl openmpi \
	--trace mpi,ucx,osrt,cuda,nvtx \
	--force-overwrite true \
	$@
