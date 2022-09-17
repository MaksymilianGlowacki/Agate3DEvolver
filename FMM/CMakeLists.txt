cmake_minimum_required(VERSION 3.21)
project(Agate3DEvolverFMM)

set(CMAKE_CXX_STANDARD 17)

add_executable(Agate3DEvolverFMM main.cpp mesh_functions_fmm.cpp mesh_functions_fmm.h ../mesh_functions.h ../mesh_functions.cpp)
target_include_directories(Agate3DEvolverFMM PRIVATE ../3rd_party_resources)

