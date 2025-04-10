set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (MU_ARCH STREQUAL "a100")
    add_compile_definitions(__GPU_VERSION)
else()
    add_compile_definitions(__CPU_VERSION)
endif()

add_library(muphys_implementation SHARED "graupel.cpp")
#include_directories(/opt/cuda/include)
#link_libraries(nvToolsExt)
#target_link_libraries(muphys_implementation muphys_core tbb)