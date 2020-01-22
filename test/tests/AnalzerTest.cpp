//
// Created by Piotr Kubala on 16/01/2020.
//

#include <catch2/catch.hpp>
#include <catch2/trompeloeil.hpp>

#include "analyzer/Analyzer.h"

#include "mocks/AnalyzerTaskMock.h"
#include "mocks/InlineAnalyzerTaskMock.h"
#include "mocks/BulkAnalyzerTaskMock.h"
#include "mocks/OstringStreamMock.h"
#include "mocks/FileUtilsMock.h"

using trompeloeil::_;

TEST_CASE("Analzer: analyze") {
    auto task1 = std::make_unique<AnalyzerTaskMock>();
    auto task2 = std::make_unique<AnalyzerTaskMock>();
    auto task3 = std::make_unique<AnalyzerTaskMock>();
    auto eigensystem = Eigensystem({1, 2, 3});
    trompeloeil::sequence seq;
    REQUIRE_CALL(*task1, analyze(eigensystem))
            .IN_SEQUENCE(seq);
    REQUIRE_CALL(*task2, analyze(eigensystem))
            .IN_SEQUENCE(seq);
    REQUIRE_CALL(*task3, analyze(eigensystem))
            .IN_SEQUENCE(seq);
    Analyzer analyzer;
    analyzer.addTask(std::move(task1));
    analyzer.addTask(std::move(task2));
    analyzer.addTask(std::move(task3));

    analyzer.analyze(eigensystem);
}

TEST_CASE("Analzyer: print inline header") {
    auto task1 = std::make_unique<InlineAnalyzerTaskMock>();
    auto task2 = std::make_unique<InlineAnalyzerTaskMock>();
    trompeloeil::sequence seq;
    REQUIRE_CALL(*task1, getResultHeader())
            .IN_SEQUENCE(seq)
            .RETURN(std::vector<std::string>{"result1", "result2"});
    REQUIRE_CALL(*task2, getResultHeader())
            .IN_SEQUENCE(seq)
            .RETURN(std::vector<std::string>{"result3"});
    Analyzer analyzer;
    analyzer.addTask(std::move(task1));
    analyzer.addTask(std::move(task2));

    auto header = analyzer.getInlineResultsHeader();

    REQUIRE(header == std::vector<std::string>{"result1", "result2", "result3"});
}

TEST_CASE("Analzyer: print inline fields") {
    auto task1 = std::make_unique<InlineAnalyzerTaskMock>();
    auto task2 = std::make_unique<InlineAnalyzerTaskMock>();
    trompeloeil::sequence seq;
    REQUIRE_CALL(*task1, getResultFields())
            .IN_SEQUENCE(seq)
            .RETURN(std::vector<std::string>{"result1", "result2"});
    REQUIRE_CALL(*task2, getResultFields())
            .IN_SEQUENCE(seq)
            .RETURN(std::vector<std::string>{"result3"});
    Analyzer analyzer;
    analyzer.addTask(std::move(task1));
    analyzer.addTask(std::move(task2));

    auto header = analyzer.getInlineResultsFields();

    REQUIRE(header == std::vector<std::string>{"result1", "result2", "result3"});
}

TEST_CASE("Analyzer: store") {
    auto out1ptr = new OstringStreamMock("task1 result");
    auto out2ptr = new OstringStreamMock("task2 result");

    auto ostreamProvider = std::make_unique<FileOstreamProviderMock>();
    REQUIRE_CALL(*ostreamProvider, openFile("sim_task1.txt"))
            .LR_RETURN(std::unique_ptr<std::ostream>(out1ptr));
    REQUIRE_CALL(*ostreamProvider, openFile("sim_task2.txt"))
            .LR_RETURN(std::unique_ptr<std::ostream>(out2ptr));

    trompeloeil::sequence seq;
    auto task1 = std::make_unique<BulkAnalyzerTaskMock>();
    REQUIRE_CALL(*task1, storeResult(_))
            .LR_WITH(&_1 == out1ptr)
            .SIDE_EFFECT(_1 << "task1 result")
            .IN_SEQUENCE(seq);
    ALLOW_CALL(*task1, getName())
            .RETURN("task1");
    auto task2 = std::make_unique<BulkAnalyzerTaskMock>();
    REQUIRE_CALL(*task2, storeResult(_))
            .LR_WITH(&_1 == out2ptr)
            .SIDE_EFFECT(_1 << "task2 result")
            .IN_SEQUENCE(seq);
    ALLOW_CALL(*task2, getName())
            .RETURN("task2");

    Analyzer analyzer(std::move(ostreamProvider));
    analyzer.addTask(std::move(task1));
    analyzer.addTask(std::move(task2));


    analyzer.storeBulkResults("sim");
}

