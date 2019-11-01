project(mbl_ed_test)

add_executable(tests test_main.cpp FockBaseGeneratorTest.cpp)
target_link_libraries(tests PRIVATE mbl_ed_src catch)