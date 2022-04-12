cmake_minimum_required(VERSION 3.21)
project(Agate3DEvolverGPU CUDA)

set(CMAKE_CUDA_STANDARD 17)

add_executable(Agate3DEvolverGPU main.cu ../mesh_functions.cpp ../mesh_functions.h)
target_include_directories(Agate3DEvolverGPU PRIVATE /Users/maks/src.git/happly)