## Implementation

### Repository
---

The project repository is located on **Levante** at the following path:

```bash
/home/b/b383366/sky/scc_at_isc25
```

### Scripts
---
`build_cpu_double.sh`: Builds the stdpar-based CPU version with double precision.

`build_cpu_single.sh`: Builds the stdpar-based CPU version with single precision.

`build_gpu_double.sh`: Builds the stdpar-based GPU+MPI version with double precision.

`build_gpu_single.sh`: Builds the stdpar-based GPU+MPI version with single precision.

`levante-gpu-energy.batch`: Specifies the optimal energy configuration (GPUs/nodes) for running the atm_R2B08.nc test cases and uses `run_wrapper_levante.sh` to generate `nvsmi.log.jobid.nodenumber`

### Logs
---

- `logs/correctness/std_<cpu|gpu>_<precision>.log`  
  Contains correctness test results for the input files `11k.nc`, `20k.nc`, and `dbg.nc`.

- `logs/optimization_strategy/gpu_<optimization_strategy>_v<version>.log`  
  Contains execution logs for each GPU optimization strategy and version, using the `/work/ka1273/atm_R2B08.nc` test case.

- `logs/energy/gpu_energy_<nodes>_<gpus_or_tasks_per_node>.log`  
  Contains energy consumption logs from multi-node executions, based on the `/work/ka1273/atm_R2B08.nc` test case.

- `logs/energy/nvsmi/nvsmi.log.<jobid>.<nodenumber>`  
  Logs collected via `nvidia-smi`, providing node-level GPU and power usage information during execution of the `/work/ka1273/atm_R2B08.nc` test case.


### Available compile options 
---

* _Implementation_
  * MU_IMPL=seq - C++ serial implementation
 * _Precision_ (default is `double`)
  * MU_ENABLE_SINGLE - to switch to `float` 
* _Enable MPI library_
    * MU_ENABLE_MPI - enable mpi (default is `OFF`)

### Modify content
---

This project integrates support for both MPI-based parallel I/O and `std::execution::par` (stdpar) parallelism. Key implementation details are outlined below:

#### 1. stdpar Implementation
- **Source file**: `<repository>/implementations/std/graupel.cpp`  
- Implements GPU-accelerated computations using C++ standard parallelism (`std::execution::par`).

#### 2. MPI Communication Type
- **File**: `<repository>/core/common/types.hpp`  
- Introduced a new type (either `float` or `double`) to support MPI-based data exchange.

#### 3. Optimized I/O for MPI Parallelism
- **Files**:  
  - `<repository>/io/io.cpp`  
  - `<repository>/io/io.hpp`  
- Added support for NetCDF parallel I/O using `NetCDF_par` to enable efficient reading and writing of files in MPI environments.

#### 4. Dual Main Functions
- **MPI-specific main function**: `<repository>/io/main_mpi.cpp`  
- When building with `MU_ENABLE_MPI=ON`, CMake selects this version to compile the `graupel` executable with MPI support.

### Usage
---
The `graupel` executable can be run in either MPI mode or serial mode depending on the build configuration.

#### MPI Mode (`MU_ENABLE_MPI=ON`)

When MPI support is enabled, use `srun` (or `mpirun`) to execute with multiple processes:

```bash
srun -n <num_procs> ./<build-dir>/bin/graupel tasks/<input-file.nc> <output-file.nc>
```

#### Serial Mode (`MU_ENABLE_MPI=OFF`) *(default)*

When MPI support is not enabled, you can directly execute the program:

```bash
./<build-dir>/bin/graupel tasks/<input-file.nc> <output-file.nc>
```


### Optimization Strategies
---

1. **Coalesced Access**  
   Optimizes GPU memory access by ensuring contiguous memory access within warps, reducing memory transactions and improving throughput. Effective in bandwidth-bound scenarios.

2. **Loop Flattening**  
   Converts nested loops into a single linear iteration space, exposing more parallelism and improving workload distribution. **Delivers the highest performance gain** among strategies by increasing occupancy and SIMD efficiency.

3. **Branchless**  
   Replaces conditional branches with arithmetic expressions (e.g., using booleans). Reduces warp divergence, improves throughput, and ensures more predictable execution, especially in tight loops.

4. **SoA (Structure of Arrays)**  
   Restructures data from AoS to separate arrays, improving memory coalescing and access patterns. Enhances GPU performance by aligning data and maximizing memory bandwidth.

5. **Pointer Capture**  
   Captures pointers directly in lambdas, avoiding unnecessary data copying. Improves memory access, reduces overhead, and enhances GPU performance by optimizing memory locality and coalescing.

### Other Optimization
---

- **IO Optimization**  
    MPI communication was found to be time-consuming, and since the data columns are independent, each rank can parallelize file reading. As a result, parallel file reading was implemented in `io.cpp` and `io.hpp` specifically for `netcdf_par`.


- **Execution Policy Optimization**

    I chose to use only `std::execution::par_unseq` for parallel execution because it provides the best performance by combining both parallel execution and SIMD (vectorization). This policy allows for efficient data processing across multiple threads while utilizing hardware vector units for data-parallel tasks, making it particularly effective for the GPU and multi-core CPUs. It offers significant performance improvements without the complexity of manually managing parallelism or vectorization.

