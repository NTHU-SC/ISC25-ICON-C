cmake_minimum_required(VERSION 3.20)
project(muphys_cpp LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
# notify all possible issues 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++20")

option(MU_ENABLE_SINGLE "Enable single precision" OFF)
# available front-ends
option(MU_IMPL "Select implementation: seq vs std" "seq")

option(MU_ENABLE_TESTS "Enable unit-tests" ON)

option(MU_ENABLE_MPI "Enable MPI support" OFF)

set(MU_ARCH "x86_64" CACHE STRING "Select architecture, x86_64, a100")

# includes
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# configure the structure of the output folders
include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}" CACHE PATH
        "Directory for the built binaries")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
        "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
        "Directory for the built libraries")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
        "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
        "Directory for the built static libraries")

if (MU_ENABLE_SINGLE)
    add_compile_definitions(__SINGLE_PRECISION)
endif ()

# add local sources
add_subdirectory(core)
add_subdirectory(io)
include(implementations/CMakeLists.txt)

# add test (if enabled)
if (MU_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(extern/googletest)
    add_subdirectory(test)
endif ()

# MPI configuration
if (MU_ENABLE_MPI)
    find_package(MPI REQUIRED)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    add_compile_definitions(USE_MPI)
    set(EXECUTABLE_SRC "main_mpi.cpp")
else()
    set(EXECUTABLE_SRC "main.cpp")
endif ()

# Build executable
add_executable(graupel ${EXECUTABLE_SRC})


# link dependency libs
target_link_libraries(graupel muphys_core muphys_io muphys_implementation)
target_include_directories(graupel PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/core"
                          "${PROJECT_SOURCE_DIR}/implementations"
                          "${PROJECT_SOURCE_DIR}/io"
                        )
set_target_properties(graupel PROPERTIES LINKER_LANGUAGE CXX)