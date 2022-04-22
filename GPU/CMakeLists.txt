cmake_minimum_required(VERSION 3.21)
project(Agate3DEvolverGPU CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(Agate3DEvolverGPU main.cu mesh_functions_gpu.cu mesh_functions_gpu.cuh)
target_include_directories(Agate3DEvolverGPU PRIVATE ../3rd_party_resources)
set_target_properties(Agate3DEvolverGPU PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES 86)
target_compile_options(Agate3DEvolverGPU PRIVATE -Xcompiler /std:c++17)