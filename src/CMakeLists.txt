project(mbl_ed_app)

add_library(mbl_ed_src STATIC Assertions.h FockBaseGenerator.cpp FockBase.cpp HamiltonianGenerator.cpp
                              CavityHamiltonianGenerator.h GapRatioCalculator.cpp Quantity.cpp Simulation.h Config.cpp
                              Utils.cpp Parameters.cpp Parameters.h)
target_include_directories(mbl_ed_src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(mbl_ed Main.cpp)
target_link_libraries(mbl_ed PRIVATE mbl_ed_src)