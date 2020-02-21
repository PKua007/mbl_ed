project(mbl_ed_test)

add_executable(tests TestMain.cpp tests/FockBaseGeneratorTest.cpp tests/FockBaseTest.cpp
        tests/HamiltonianGeneratorTest.cpp tests/MeanGapRatioTest.cpp tests/SimulationTest.cpp
        matchers/ArmaApproxEqualMatcher.h tests/AnalzerTest.cpp mocks/AnalyzerTaskMock.h  mocks/FileUtilsMock.h
        mocks/OstringStreamMock.h mocks/InlineAnalyzerTaskMock.h mocks/BulkAnalyzerTaskMock.h tests/CDFTest.cpp
        tests/FoldTest.cpp tests/EigensystemTest.cpp tests/MeanInverseParticipationRatioTest.cpp
        tests/InverseParticipationRatioTest.cpp tests/AveragingModelsTest.cpp mocks/RNDMock.h
        mocks/DiagonalTermMock.h mocks/HoppingTermMock.h tests/HubbardOnsite.cpp tests/CavityLongInteractionTest.cpp
        tests/OnsiteDisorderTest.cpp tests/HubbardHopTest.cpp tests/ConfigTest.cpp tests/CavityConstantsTest.cpp
        tests/CavityConstantsReaderTest.cpp tests/LookupCavityZ2Test.cpp tests/LookupCavityYZTest.cpp
        mocks/DoubleHoppingTermMock.h tests/LookupCavityY2Test.cpp tests/CorrelationsTimeEvolutionTest.cpp
        tests/SymmetricMatrixTest.cpp)
target_link_libraries(tests PRIVATE mbl_ed_src Catch2::Catch2 trompeloeil)
target_include_directories(tests PRIVATE ../test)