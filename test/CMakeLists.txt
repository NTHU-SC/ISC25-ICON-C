include(FetchContent)

enable_testing()

# include all source files
include_directories(${CMAKE_SOURCE_DIR})
   
# create test executable
add_executable(muphys_core_test common.cc)
if (MU_ENABLE_STANDALONE)
  target_sources(muphys_core_test PRIVATE io.cc)
endif()

target_include_directories(muphys_core_test PRIVATE ${CMAKE_SOURCE_DIR})

# link against googletest (built locally from /extern)
target_link_libraries(muphys_core_test GTest::gtest_main muphys_core muphys_io muphys_implementation)

# LLVM requires filesystem library linked explicitly
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "NVHPC")
  target_link_libraries(muphys_core_test stdc++fs)
endif()  

include(GoogleTest)
gtest_discover_tests(muphys_core_test)