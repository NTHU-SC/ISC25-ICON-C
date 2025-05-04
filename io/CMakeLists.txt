# find netcdf
find_package(NetCDF COMPONENTS C CXX REQUIRED)
include_directories(${NetCDF_CXX_INCLUDE_DIR})

if(MU_ENABLE_MPI)
find_package(MPI REQUIRED)
add_compile_definitions(USE_MPI)
include_directories(${MPI_INCLUDE_PATH})
endif()

add_library(muphys_io SHARED "io.cpp")
set_target_properties(muphys_io PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(muphys_io PUBLIC NetCDF::NetCDF_CXX NetCDF::NetCDF_C)

if(MU_ENABLE_MPI)
target_include_directories(muphys_io PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(muphys_io PRIVATE ${MPI_CXX_LIBRARIES})
endif()