#! /bin/bash
rm -rf nsys_reports
mkdir -p nsys_reports

# Output to ./nsys_reports/rank_$N.nsys-rep
/opt/nvidia/nsight-systems/2025.1.1/bin/nsys  profile \
	-o "./nsys_reports/gpu_${OMPI_COMM_WORLD_LOCAL_RANK}.nsys-rep" \
	--mpi-impl openmpi \
	--trace mpi,ucx,osrt,cuda,nvtx \
	--force-overwrite true \
	$@
