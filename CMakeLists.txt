cmake_minimum_required(VERSION 3.15)
project(cavity_mbl)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

add_executable(cavity_mbl main.cpp)
target_link_libraries(cavity_mbl -larmadillo)