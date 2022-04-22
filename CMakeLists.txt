cmake_minimum_required(VERSION 3.21)
project(Agate3DEvolver)

set(CMAKE_CXX_STANDARD 17)

if(CMAKE_BUILD_TYPE STREQUAL Release)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/")
endif()

add_executable(Agate3DEvolver main.cpp mesh_functions.cpp mesh_functions.h)
target_include_directories(Agate3DEvolver PRIVATE 3rd_party_resources)

if(WIN32)
    add_subdirectory(GPU)
    add_subdirectory(FIM)
    add_subdirectory(FMM)
endif()