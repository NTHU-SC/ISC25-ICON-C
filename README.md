## Installation

### Dependencies
* [NetCDF for CXX](https://github.com/Unidata/netcdf-cxx4)
  * for Levante: `spack load netcdf-cxx4@4.3.1`
* nvc++ compiler
  * for Levante `module load nvhpc/24.7-gcc-11.2.0`
* MPI library (& compiler)
  * for Levante `module load openmpi/4.1.5-nvhpc-24.7` and the compiler to use is then `mpicxx`
For automatic setup, use `source scripts/levante-setup.sh`.


### Available compile options 
* _Implementation_ - The sequential implementation is selected by default. The user can choose of the following options:
  * MU_IMPL=seq - C++ serial implementation
 * _Precision_ (default is `double`)
  * MU_ENABLE_SINGLE - to switch to `float` 
* _Unit-test_ - compile tests together with the main executable (default is `true`)
  * MU_ENABLE_TESTS
  * MU_ENABLE_MPI - enable mpi (default is `OFF`)

### Compile the project (with default flags and Seq frontend)

`cmake -DMU_IMPL=seq -B <build-dir> -S . -DCMAKE_CXX_COMPILER=nvc++`

`cmake --build <build-dir>`

## Usage

`./<build-dir>/bin/graupel tasks/<input-file.nc> <output-file.nc>`

## Automated tests

- Run tests manually:
`cd <build-dir> && ctest` 

