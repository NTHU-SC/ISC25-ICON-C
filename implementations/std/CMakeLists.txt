add_library(muphys_implementation SHARED "graupel.cpp")
include_directories(/opt/cuda/include)
link_libraries(nvToolsExt)
target_link_libraries(muphys_implementation muphys_core tbb)