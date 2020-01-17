//
// Created by Piotr Kubala on 16/01/2020.
//

#include <catch2/catch.hpp>
#include <catch2/trompeloeil.hpp>

#include "Analyzer.h"

#include "AnalyzerTaskMock.h"
#include "OstringStreamMock.h"
#include "FileUtilsMock.h"

using trompeloeil::_;

TEST_CASE("Analzer: analyze") {
    auto task1 = std::make_unique<AnalyzerTaskMock>();
    auto task2 = std::make_unique<AnalyzerTaskMock>();
    auto task3 = std::make_unique<AnalyzerTaskMock>();
    std::vector<double> eigenenergies = {1, 2, 3};
    trompeloeil::sequence seq;
    REQUIRE_CALL(*task1, analyze(eigenenergies))
            .IN_SEQUENCE(seq);
    REQUIRE_CALL(*task2, analyze(eigenenergies))
            .IN_SEQUENCE(seq);
    REQUIRE_CALL(*task3, analyze(eigenenergies))
            .IN_SEQUENCE(seq);
    Analyzer analyzer;
    analyzer.addTask(std::move(task1));
    analyzer.addTask(std::move(task2));
    analyzer.addTask(std::move(task3));

    analyzer.analyze(eigenenergies);
}

TEST_CASE("Analzer: print") {
    auto task1 = std::make_unique<AnalyzerTaskMock>();
    auto task2 = std::make_unique<AnalyzerTaskMock>();
    std::ostringstream out;
    trompeloeil::sequence seq;
    using trompeloeil::_;
    REQUIRE_CALL(*task1, printResult(_))
            .LR_WITH(&_1 == &out)
            .SIDE_EFFECT(_1 << "task1 result" << std::endl)
            .IN_SEQUENCE(seq);
    REQUIRE_CALL(*task2, printResult(_))
            .LR_WITH(&_1 == &out)
            .SIDE_EFFECT(_1 << "task2 result" << std::endl)
            .IN_SEQUENCE(seq);
    REQUIRE_CALL(*task1, getName())
            .RETURN("task1");
    REQUIRE_CALL(*task2, getName())
            .RETURN("task2");
    Analyzer analyzer;
    analyzer.addTask(std::move(task1));
    analyzer.addTask(std::move(task2));

    analyzer.printResults(out);

    std::ostringstream expected;
    expected << "task1 results:" << std::endl << "task1 result" << std::endl;
    expected << "task2 results:" << std::endl << "task2 result" << std::endl;
    REQUIRE(out.str() == expected.str());
}

TEST_CASE("Analzer: store") {
    auto out1ptr = new OstringStreamMock("task1 result\n");
    auto out2ptr = new OstringStreamMock("task2 result\n");

    auto ostreamProvider = std::make_unique<FileOstreamProviderMock>();
    REQUIRE_CALL(*ostreamProvider, openFile("sim_task1.txt"))
            .LR_RETURN(std::unique_ptr<std::ostream>(out1ptr));
    REQUIRE_CALL(*ostreamProvider, openFile("sim_task2.txt"))
            .LR_RETURN(std::unique_ptr<std::ostream>(out2ptr));

    trompeloeil::sequence seq;
    auto task1 = std::make_unique<AnalyzerTaskMock>();
    REQUIRE_CALL(*task1, printResult(_))
            .LR_WITH(&_1 == out1ptr)
            .SIDE_EFFECT(_1 << "task1 result" << std::endl)
            .IN_SEQUENCE(seq);
    ALLOW_CALL(*task1, getName())
            .RETURN("task1");
    auto task2 = std::make_unique<AnalyzerTaskMock>();
    REQUIRE_CALL(*task2, printResult(_))
            .LR_WITH(&_1 == out2ptr)
            .SIDE_EFFECT(_1 << "task2 result" << std::endl)
            .IN_SEQUENCE(seq);
    ALLOW_CALL(*task2, getName())
            .RETURN("task2");

    Analyzer analyzer(std::move(ostreamProvider));
    analyzer.addTask(std::move(task1));
    analyzer.addTask(std::move(task2));


    analyzer.storeResults("sim");
}

