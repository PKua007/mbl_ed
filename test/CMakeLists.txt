project(mbl_ed_test)

add_executable(tests TestMain.cpp tests/FockBaseGeneratorTest.cpp tests/FockBaseTest.cpp tests/HamiltonianGeneratorTest.cpp
        tests/CavityHamiltonianGeneratorTest.cpp tests/GapRatioCalculatorTest.cpp tests/SimulationTest.cpp
        matchers/ArmaMatrixApproxEqualMatcher.h tests/AnalzerTest.cpp mocks/AnalyzerTaskMock.h mocks/FileUtilsMock.h
        mocks/OstringStreamMock.h mocks/InlineAnalyzerTaskMock.h mocks/BulkAnalyzerTaskMock.h)
target_link_libraries(tests PRIVATE mbl_ed_src Catch2::Catch2 trompeloeil)
target_include_directories(tests PRIVATE ../test)