cmake_minimum_required(VERSION 3.21)
project(Agate3DEvolverFMM)

set(CMAKE_CXX_STANDARD 17)

add_executable(Agate3DEvolverFMM main.cpp)
target_include_directories(Agate3DEvolverFMM PRIVATE ../3rd_party_resources)

