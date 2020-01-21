project(mbl_ed_app)

add_library(mbl_ed_src STATIC utils/Assertions.h simulation/FockBaseGenerator.cpp simulation/FockBase.cpp simulation/HamiltonianGenerator.cpp
        simulation/CavityHamiltonianGenerator.h analyzer/tasks/MeanGapRatio.cpp utils/Quantity.cpp simulation/Simulation.h utils/Config.cpp
        utils/Utils.cpp Parameters.cpp Parameters.h analyzer/AnalyzerTask.h analyzer/Analyzer.cpp analyzer/Analyzer.h utils/FileUtils.h analyzer/InlineAnalyzerTask.h analyzer/BulkAnalyzerTask.h analyzer/tasks/CDF.cpp analyzer/tasks/CDF.h InlineResultsPrinter.cpp InlineResultsPrinter.h utils/Fold.cpp utils/Fold.h)
target_include_directories(mbl_ed_src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(mbl_ed Main.cpp)
target_link_libraries(mbl_ed PRIVATE mbl_ed_src cxxopts)