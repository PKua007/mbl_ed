project(mbl_ed_app)

add_library(mbl_ed_src STATIC utils/Assertions.h core/FockBasisGenerator.cpp core/FockBasis.cpp
        core/HamiltonianGenerator.cpp analyzer/tasks/MeanGapRatio.cpp
        utils/Quantity.cpp simulation/ExactDiagonalization.h utils/Config.cpp utils/Utils.cpp frontend/Parameters.cpp
        analyzer/AnalyzerTask.h analyzer/Analyzer.cpp utils/FileUtils.h analyzer/InlineAnalyzerTask.h
        analyzer/BulkAnalyzerTask.h analyzer/tasks/CDF.cpp utils/Fold.cpp
        frontend/IO.cpp frontend/Frontend.cpp core/Eigensystem.cpp
        analyzer/tasks/MeanInverseParticipationRatio.cpp simulation/ExactDiagonalizationParameters.h
        analyzer/tasks/InverseParticipationRatio.cpp analyzer/tasks/InverseParticipationRatio.h
        core/AveragingModel.h core/DisorderGenerator.h core/RND.h core/DiagonalTerm.h
        core/HoppingTerm.h core/terms/HubbardHop.cpp core/terms/HubbardOnsite.cpp
        core/terms/OnsiteDisorder.h core/terms/CavityLongInteraction.cpp core/CavityConstants.cpp
        frontend/CavityConstantsReader.cpp core/terms/LookupCavityZ2.cpp core/terms/LookupCavityYZ.cpp
        core/DoubleHoppingTerm.h core/terms/LookupCavityY2.cpp
        analyzer/tasks/EDTimeEvolution.cpp evolution/OservablesTimeEvolution.cpp evolution/SymmetricMatrix.h
        evolution/TimeEvolutionEntry.cpp core/terms/QuasiperiodicDisorder.cpp
        core/terms/QuasiperiodicDisorder.h core/terms/ListOnsite.cpp evolution/Evolver.h
        evolution/EDEvolver.cpp evolution/ChebyshevEvolver.cpp evolution/TimeEvolution.cpp
        simulation/ChebyshevEvolution.h evolution/TimeEvolutionParameters.h
        evolution/EvolutionTimeSegment.h frontend/HamiltonianGeneratorBuilder.cpp frontend/AnalyzerBuilder.cpp
        core/terms/OnsiteDisorder.cpp frontend/AveragingModelFactory.cpp
        core/averaging_models/DummyAveragingModel.h core/averaging_models/CavityConstantsAveragingModel.cpp
        core/averaging_models/OnsiteDisorderAveragingModel.cpp
        core/averaging_models/RandomPhi0AveragingModel.cpp
        core/averaging_models/UniformPhi0AveragingModel.cpp
        core/disorder_generators/UniformGenerator.h core/FockVector.cpp
        analyzer/tasks/DressedStatesFinder.cpp analyzer/tasks/BulkMeanGapRatio.cpp core/QuenchCalculator.cpp
        simulation/SimulationsSpan.h utils/OMPMacros.h simulation/QuenchDataSimulation.h simulation/Restorable.h
        simulation/RestorableSimulation.h simulation/RestorableSimulationExecutor.cpp simulation/RestorableHelper.h
        utils/Logger.h core/Observable.h core/PrimaryObservable.h core/SecondaryObservable.h
        evolution/TimeEvolutionParameters.cpp core/observables/OnsiteOccupations.cpp
        core/observables/OnsiteFluctuations.cpp core/observables/Correlations.cpp
        core/observables/OnsiteOccupationsSquared.cpp frontend/ObservablesBuilder.cpp
        core/observables/BipariteEntropy.cpp analyzer/tasks/EigenstateObservables.cpp analyzer/tasks/PDF.cpp
        simulation/RandomStateObservables.cpp core/observables/CavityOnsiteOccupations.cpp
        core/observables/CavityOnsiteOccupationsSquared.cpp core/observables/CavityElectricField.cpp
        core/observables/CavityLightIntensity.cpp analyzer/tasks/ParticipationEntropy.cpp
        analyzer/EigenvaluesExtractor.cpp)

target_include_directories(mbl_ed_src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(mbl_ed_src PUBLIC ../extern/ZipIterator)
target_link_libraries(mbl_ed_src PUBLIC cxxopts stdc++fs)

add_executable(mbl_ed Main.cpp)
target_link_libraries(mbl_ed PRIVATE mbl_ed_src)