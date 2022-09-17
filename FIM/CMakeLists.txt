cmake_minimum_required(VERSION 3.21)
project(Agate3DEvolverFIM CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(Agate3DEvolverFIM main.cu ../FMM/mesh_functions_fmm.cpp ../FMM/mesh_functions_fmm.h ../mesh_functions.h ../mesh_functions.cpp)
target_link_libraries(Agate3DEvolverFIM "C:/Programming/SCI_Solver/StructuredEikonal/build/STRUCTURED_EIKONAL.lib")
target_include_directories(Agate3DEvolverFIM PRIVATE "C:/Programming/SCI_Solver/StructuredEikonal/src" ../3rd_party_resources)

target_compile_options(Agate3DEvolverFIM PRIVATE -Xcompiler /std:c++17)
set_target_properties(Agate3DEvolverFIM PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES 86)
