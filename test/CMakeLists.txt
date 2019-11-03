project(mbl_ed_test)

add_executable(tests TestMain.cpp FockBaseGeneratorTest.cpp FockBaseTest.cpp HamiltonianGeneratorTest.cpp)
target_link_libraries(tests PRIVATE mbl_ed_src catch)