project(mbl_ed_app)

add_library(mbl_ed_src STATIC FockBaseGenerator.cpp)
target_include_directories(mbl_ed_src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mbl_ed_src INTERFACE -larmadillo)

add_executable(mbl_ed main.cpp)
target_link_libraries(mbl_ed PRIVATE mbl_ed_src)