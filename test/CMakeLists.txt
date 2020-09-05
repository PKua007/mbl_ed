project(mbl_ed_test)

add_executable(tests TestMain.cpp tests/core/FockBaseGeneratorTest.cpp tests/core/FockBaseTest.cpp
        tests/core/HamiltonianGeneratorTest.cpp tests/analyzer/MeanGapRatioTest.cpp
        tests/simulation/ExactDiagonalizationTest.cpp matchers/ArmaApproxEqualMatcher.h
        tests/analyzer/AnalzerTest.cpp mocks/AnalyzerTaskMock.h  mocks/FileUtilsMock.h mocks/OstringStreamMock.h
        mocks/InlineAnalyzerTaskMock.h mocks/BulkAnalyzerTaskMock.h tests/analyzer/CDFTest.cpp
        tests/utils/FoldTest.cpp tests/core/EigensystemTest.cpp tests/analyzer/MeanInverseParticipationRatioTest.cpp
        tests/analyzer/InverseParticipationRatioTest.cpp tests/core/AveragingModelsTest.cpp mocks/RNDMock.h
        mocks/DiagonalTermMock.h mocks/HoppingTermMock.h tests/core/HubbardOnsite.cpp
        tests/core/CavityLongInteractionTest.cpp tests/core/OnsiteDisorderTest.cpp tests/core/HubbardHopTest.cpp
        tests/utils/ConfigTest.cpp tests/core/CavityConstantsTest.cpp tests/frontend/CavityConstantsReaderTest.cpp
        tests/core/LookupCavityZ2Test.cpp tests/core/LookupCavityYZTest.cpp mocks/DoubleHoppingTermMock.h
        tests/core/LookupCavityY2Test.cpp tests/evolution/CorrelationsTimeEvolutionTest.cpp
        tests/evolution/SymmetricMatrixTest.cpp tests/evolution/OccupationEvolutionTest.cpp
        matchers/VectorApproxEqualMatcher.h  tests/evolution/CorrelationsTimeEntryTest.cpp
        tests/core/QuasiperiodicDisorderTest.cpp tests/core/ListOnsiteTest.cpp tests/evolution/EvolversTest.cpp
        tests/core/FockVectorTest.cpp tests/analyzer/DressedStatesFinderTest.cpp tests/analyzer/BulkMeanGapRatioTest.cpp
        tests/core/QuenchCalculatorTest.cpp tests/simulation/ChebyshevEvolutionTest.cpp mocks/RestorableSimulationMock.h
        tests/simulation/RestorableSimulationExecutorTest.cpp object_mothers/HamiltonianGeneratorMother.cpp tests/simulation/RestorableContractTest.cpp)
target_link_libraries(tests PRIVATE mbl_ed_src Catch2::Catch2 trompeloeil)
target_include_directories(tests PRIVATE ../test)