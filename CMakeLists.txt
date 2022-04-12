cmake_minimum_required(VERSION 3.21)
project(Agate3DEvolver)

set(CMAKE_CXX_STANDARD 17)

add_executable(Agate3DEvolver main.cpp mesh_functions.cpp mesh_functions.h)
target_include_directories(Agate3DEvolver PRIVATE /Users/maks/src.git/happly)
