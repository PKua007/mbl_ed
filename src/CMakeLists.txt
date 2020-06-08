project(mbl_ed_app)

add_library(mbl_ed_src STATIC utils/Assertions.h simulation/FockBaseGenerator.cpp simulation/FockBase.cpp
        simulation/HamiltonianGenerator.cpp analyzer/tasks/MeanGapRatio.cpp
        utils/Quantity.cpp simulation/Simulation.h utils/Config.cpp utils/Utils.cpp frontend/Parameters.cpp
        analyzer/AnalyzerTask.h analyzer/Analyzer.cpp utils/FileUtils.h analyzer/InlineAnalyzerTask.h
        analyzer/BulkAnalyzerTask.h analyzer/tasks/CDF.cpp utils/Fold.cpp
        frontend/IO.cpp frontend/Frontend.cpp simulation/Eigensystem.cpp
        analyzer/tasks/MeanInverseParticipationRatio.cpp simulation/SimulationParameters.h
        analyzer/tasks/InverseParticipationRatio.cpp analyzer/tasks/InverseParticipationRatio.h
        simulation/AveragingModel.h simulation/DisorderGenerator.h simulation/RND.h simulation/DiagonalTerm.h
        simulation/HoppingTerm.h simulation/terms/HubbardHop.cpp simulation/terms/HubbardOnsite.cpp
        simulation/terms/OnsiteDisorder.h simulation/terms/CavityLongInteraction.cpp simulation/CavityConstants.cpp
        frontend/CavityConstantsReader.cpp simulation/terms/LookupCavityZ2.cpp simulation/terms/LookupCavityYZ.cpp
        simulation/DoubleHoppingTerm.h simulation/terms/LookupCavityY2.cpp
        analyzer/tasks/EDCorrelationsTimeEvolution.cpp evolution/OccupationEvolution.cpp evolution/SymmetricMatrix.h
        evolution/CorrelationsTimeEntry.cpp simulation/terms/QuasiperiodicDisorder.cpp
        simulation/terms/QuasiperiodicDisorder.h simulation/terms/ListOnsite.cpp evolution/Evolver.h
        evolution/EDEvolver.cpp evolution/ChebyshevEvolver.cpp evolution/CorrelationsTimeEvolution.cpp
        evolution/ChebyshevEvolution.h evolution/CorrelationsTimeEvolutionParameters.h evolution/EvolutionTimeSegment.h frontend/HamiltonianGeneratorBuilder.cpp frontend/HamiltonianGeneratorBuilder.h frontend/AnalyzerBuilder.cpp frontend/AnalyzerBuilder.h simulation/terms/OnsiteDisorder.cpp frontend/AveragingModelFactory.cpp frontend/AveragingModelFactory.h)
target_include_directories(mbl_ed_src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(mbl_ed_src PUBLIC ../extern/ZipIterator)
target_link_libraries(mbl_ed_src PUBLIC cxxopts stdc++fs)

add_executable(mbl_ed Main.cpp simulation/averaging_models/DummyAveragingModel.h simulation/averaging_models/OnsiteDisorderAveragingModel.h simulation/averaging_models/UniformPhi0AveragingModel.h simulation/averaging_models/RandomPhi0AveragingModel.h simulation/averaging_models/CavityConstantsAveragingModel.h simulation/averaging_models/CavityConstantsAveragingModel.cpp simulation/averaging_models/OnsiteDisorderAveragingModel.cpp simulation/averaging_models/RandomPhi0AveragingModel.cpp simulation/averaging_models/UniformPhi0AveragingModel.cpp simulation/disorder_generators/UniformGenerator.h)
target_link_libraries(mbl_ed PRIVATE mbl_ed_src)