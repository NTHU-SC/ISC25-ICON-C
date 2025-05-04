## Implementation

### Repository
---
```
/home/b/b383366/sky/scc_at_isc25
```

### Build Script
---
`build_cpu_double.sh`: Builds the stdpar-based CPU version with double precision.

`build_cpu_single.sh`: Builds the stdpar-based CPU version with single precision.

`build_gpu_double.sh`: Builds the stdpar-based GPU+MPI version with double precision.

`build_gpu_single.sh`: Builds the stdpar-based GPU+MPI version with single precision.

### Dependencies
---

* [NetCDF for CXX](https://github.com/Unidata/netcdf-cxx4)
  * for Levante: `spack load netcdf-cxx4@4.3.1`
* nvc++ compiler
  * for Levante `module load nvhpc/24.7-gcc-11.2.0`
* MPI library (& compiler)
  * for Levante `module load openmpi/4.1.5-nvhpc-24.7` and the compiler to use is then `mpicxx`
For automatic setup, use `source scripts/levante-setup.sh`.


### Available compile options 
---

* _Implementation_ - The sequential implementation is selected by default. The user can choose of the following options:
  * MU_IMPL=seq - C++ serial implementation
 * _Precision_ (default is `double`)
  * MU_ENABLE_SINGLE - to switch to `float` 
* _Unit-test_ - compile tests together with the main executable (default is `true`)
  * MU_ENABLE_TESTS
  * MU_ENABLE_MPI - enable mpi (default is `OFF`)
