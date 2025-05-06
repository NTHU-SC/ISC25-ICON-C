# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/b/b383366/sky/scc_at_isc25/build_std_cpu_double/_deps/googletest-src"
  "/home/b/b383366/sky/scc_at_isc25/build_std_cpu_double/_deps/googletest-build"
  "/home/b/b383366/sky/scc_at_isc25/build_std_cpu_double/_deps/googletest-subbuild/googletest-populate-prefix"
  "/home/b/b383366/sky/scc_at_isc25/build_std_cpu_double/_deps/googletest-subbuild/googletest-populate-prefix/tmp"
  "/home/b/b383366/sky/scc_at_isc25/build_std_cpu_double/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
  "/home/b/b383366/sky/scc_at_isc25/build_std_cpu_double/_deps/googletest-subbuild/googletest-populate-prefix/src"
  "/home/b/b383366/sky/scc_at_isc25/build_std_cpu_double/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/b/b383366/sky/scc_at_isc25/build_std_cpu_double/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp/${subDir}")
endforeach()
