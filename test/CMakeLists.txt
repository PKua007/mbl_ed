project(mbl_ed_test)

add_executable(tests TestMain.cpp FockBaseGeneratorTest.cpp FockBaseTest.cpp HamiltonianGeneratorTest.cpp
                     CavityHamiltonianGeneratorTest.cpp GapRatioCalculatorTest.cpp SimulationTest.cpp ArmaMatrixApproxEqualMatcher.h AnalzerTest.cpp AnalyzerTaskMock.h FileUtilsMock.h OstringStreamMock.h)
target_link_libraries(tests PRIVATE mbl_ed_src Catch2::Catch2 trompeloeil)