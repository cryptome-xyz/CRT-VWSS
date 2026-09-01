cmake_minimum_required(VERSION 3.16)
project(CRT_VSS LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_executable(test
    src/benchmark.cpp
    src/hash.cpp
    src/intcom.cpp
    src/vwss.cpp
    src/three_square.cpp
)

target_include_directories(test
    PRIVATE
        include
)

target_link_libraries(test
    PRIVATE
        ntl
        pari
        ssl
        crypto
)
