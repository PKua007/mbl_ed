cmake_minimum_required(VERSION 3.10)
project(mbl_ed)

set(CMAKE_CXX_STANDARD 17)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
link_libraries(${ARMADILLO_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} --coverage")

add_subdirectory(extern/Catch2)
add_subdirectory(extern/trompeloeil)
add_subdirectory(extern/cxxopts)
add_subdirectory(src)
add_subdirectory(test)