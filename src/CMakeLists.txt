project(mbl_ed_app)

add_library(mbl_ed_src STATIC utils/Assertions.h simulation/FockBaseGenerator.cpp simulation/FockBase.cpp
        simulation/HamiltonianGenerator.cpp simulation/CavityHamiltonianGenerator.h analyzer/tasks/MeanGapRatio.cpp
        utils/Quantity.cpp simulation/Simulation.h utils/Config.cpp utils/Utils.cpp frontend/Parameters.cpp
        analyzer/AnalyzerTask.h analyzer/Analyzer.cpp utils/FileUtils.h analyzer/InlineAnalyzerTask.h
        analyzer/BulkAnalyzerTask.h analyzer/tasks/CDF.cpp frontend/InlineResultsPrinter.cpp utils/Fold.cpp
        frontend/IO.cpp frontend/Frontend.cpp)
target_include_directories(mbl_ed_src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mbl_ed_src PUBLIC cxxopts stdc++fs)

add_executable(mbl_ed Main.cpp)
target_link_libraries(mbl_ed PRIVATE mbl_ed_src)