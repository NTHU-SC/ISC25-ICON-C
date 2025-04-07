export NETCDF_ROOT=/home/lin/spack/opt/spack/linux-x86_64_v4/netcdf-cxx4-4.3.1-zup3vfiehzolwgeedwprsfnw57npjyjn
export NetCDF_CXX_LIBRARY=${NETCDF_ROOT}/lib/libnetcdf-cxx4.so

cmake -DMU_IMPL=seq -DNetCDF_CXX_LIBRARY=${NetCDF_CXX_LIBRARY} -B ./build -S . -DCMAKE_CXX_COMPILER=nvc++
cmake -DMU_IMPL=seq -B ./build -S . -DCMAKE_CXX_COMPILER=nvc++
cmake --build ./build

cmake -DMU_IMPL=seq -B ./build -S . -DCMAKE_CXX_COMPILER=nvc++ \
      -DNetCDF_CXX_LIBRARY=/home/lin/spack/opt/spack/linux-x86_64_v4/netcdf-cxx4-4.3.1-zup3vfiehzolwgeedwprsfnw57npjyjn/lib/libnetcdf_c++.a \
      -DNetCDF_INCLUDE_DIR=/home/lin/spack/opt/spack/linux-x86_64_v4/netcdf-cxx4-4.3.1-zup3vfiehzolwgeedwprsfnw57npjyjn/include
