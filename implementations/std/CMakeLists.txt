set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(muphys_implementation SHARED "graupel.cpp")

if (MU_ARCH STREQUAL "x86_64")
    target_compile_options(muphys_implementation PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,NVIDIA>:--stdpar=multicore>)
    target_compile_definitions(muphys_implementation PRIVATE __CPU_VERSION)
else()
    target_compile_options(muphys_implementation PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,NVIDIA>:--stdpar>)
    target_compile_definitions(muphys_implementation PRIVATE __GPU_VERSION)
endif()

target_link_libraries(muphys_implementation PRIVATE muphys_core tbb)

#include_directories(/opt/cuda/include)
#link_libraries(nvToolsExt)
#target_link_libraries(muphys_implementation muphys_core tbb)