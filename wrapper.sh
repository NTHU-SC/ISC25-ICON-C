#! /bin/bash

# Output to ./nsys_reports/rank_$N.nsys-rep
nsys  profile \
	-o "./nsys_reports/gpu_${OMPI_COMM_WORLD_LOCAL_RANK}.nsys-rep" \
	--mpi-impl openmpi \
	--trace mpi,ucx,osrt,cuda,nvtx \
	--force-overwrite true \
	$@
