//
// Created by Piotr Kubala on 16/01/2020.
//

#include <numeric>

#include <catch2/catch.hpp>
#include <catch2/trompeloeil.hpp>

#include "analyzer/Analyzer.h"

#include "mocks/AnalyzerTaskMock.h"
#include "mocks/InlineAnalyzerTaskMock.h"
#include "mocks/BulkAnalyzerTaskMock.h"
#include "mocks/OstringStreamMock.h"
#include "mocks/FileUtilsMock.h"

#include "utils/Assertions.h"

using trompeloeil::_;

namespace {
    class MeanEnergyCalculator : public InlineAnalyzerTask {
    private:
        bool normalized{};
        std::vector<double> energies;

    public:
        explicit MeanEnergyCalculator(bool normalized) : normalized{normalized}
        { }

        void analyze(const Eigensystem &eigensystem, std::ostream &) override {
            arma::vec newEnergies;
            if (this->normalized)
                newEnergies = eigensystem.getNormalizedEigenenergies();
            else
                newEnergies = eigensystem.getEigenenergies();
            this->energies.insert(this->energies.end(), newEnergies.begin(), newEnergies.end());
        }

        std::string getName() const override { return "mean_nrg"; }
        std::vector<std::string> getResultHeader() const override { return {"mean_energy"}; }

        std::vector<std::string> getResultFields() const override {
            double mean = std::accumulate(this->energies.begin(), this->energies.end(), 0.0) / this->energies.size();
            return {std::to_string(mean)};
        }

        void storeState(std::ostream &binaryOut) const override {
            std::size_t energiesSize = this->energies.size();
            binaryOut.write(reinterpret_cast<const char*>(&energiesSize), sizeof(energiesSize));
            binaryOut.write(reinterpret_cast<const char*>(this->energies.data()),
                            sizeof(this->energies[0]) * energiesSize);
            Assert(binaryOut.good());
        }

        void joinRestoredState(std::istream &binaryIn) override {
            std::size_t energiesSizeRestored{};
            binaryIn.read(reinterpret_cast<char*>(&energiesSizeRestored), sizeof(energiesSizeRestored));
            Assert(binaryIn.good());
            std::vector<double> energiesRestored(energiesSizeRestored);
            binaryIn.read(reinterpret_cast<char*>(energiesRestored.data()),
                          sizeof(energiesRestored[0]) * energiesSizeRestored);
            Assert(binaryIn.good());
            this->energies.insert(this->energies.end(), energiesRestored.begin(), energiesRestored.end());
        }

        void clear() override {
            this->energies.clear();
        }
    };
}

TEST_CASE("Analzer: analyze") {
    auto task1 = std::make_unique<AnalyzerTaskMock>();
    auto task2 = std::make_unique<AnalyzerTaskMock>();
    auto task3 = std::make_unique<AnalyzerTaskMock>();
    auto eigensystem = Eigensystem({1, 2, 3});
    trompeloeil::sequence seq;
    REQUIRE_CALL(*task1, analyze(eigensystem, _))
            .IN_SEQUENCE(seq);
    REQUIRE_CALL(*task2, analyze(eigensystem, _))
            .IN_SEQUENCE(seq);
    REQUIRE_CALL(*task3, analyze(eigensystem, _))
            .IN_SEQUENCE(seq);
    Analyzer analyzer;
    analyzer.addTask(std::move(task1));
    analyzer.addTask(std::move(task2));
    analyzer.addTask(std::move(task3));
    std::ostringstream logger;

    analyzer.analyze(eigensystem, logger);
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
    REQUIRE_CALL(*ostreamProvider, openOutputFile("sim_task1.txt", false))
            .LR_RETURN(std::unique_ptr<std::ostream>(out1ptr));
    REQUIRE_CALL(*ostreamProvider, openOutputFile("sim_task2.txt", false))
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

TEST_CASE("Analyzer: storing, restoring and cleaning") {
    Eigensystem eigensystem0({10, 10.1, 10.15, 10.16, 10.2, 10.25, 10.47, 10.58, 10.71, 10.76, 10.8, 10.94, 11.0});
    Eigensystem eigensystem1({20, 20.05, 20.25, 20.28, 20.32, 20.58, 20.59, 20.6, 20.75, 20.85, 20.91, 21.0});
    Analyzer analyzer;
    analyzer.addTask(std::make_unique<MeanEnergyCalculator>(true));
    analyzer.addTask(std::make_unique<MeanEnergyCalculator>(false));
    std::ostringstream logger;

    SECTION("clearing") {
        analyzer.analyze(eigensystem0, logger);
        auto result1 = analyzer.getInlineResultsFields();

        analyzer.analyze(eigensystem1, logger);
        analyzer.clear();
        analyzer.analyze(eigensystem0, logger);
        auto result2 = analyzer.getInlineResultsFields();

        REQUIRE(result1 == result2);
    }

    SECTION("storing and joining restored") {
        analyzer.analyze(eigensystem0, logger);
        analyzer.analyze(eigensystem1, logger);
        auto normalResult = analyzer.getInlineResultsFields();

        analyzer.clear();
        analyzer.analyze(eigensystem1, logger);
        std::stringstream simulation1;
        analyzer.storeState(simulation1);

        analyzer.clear();
        analyzer.analyze(eigensystem0, logger);
        analyzer.joinRestoredState(simulation1);
        auto restoredResult = analyzer.getInlineResultsFields();

        REQUIRE(normalResult == restoredResult);
    }
}
