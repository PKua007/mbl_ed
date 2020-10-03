//
// Created by Piotr Kubala on 03/10/2020.
//

#include <catch2/catch.hpp>
#include <catch2/trompeloeil.hpp>

#include "matchers/TrompeloeilMatchers.h"

#include "mocks/PrimaryObservableMock.h"
#include "mocks/SecondaryObservableMock.h"

#include "analyzer/tasks/EigenstateObservables.h"

using trompeloeil::_;

namespace {
    inline auto first_observable_values_eq(const std::vector<double> &values) {
        return trompeloeil::make_matcher<std::vector<std::shared_ptr<PrimaryObservable>>>(
                [](const std::vector<std::shared_ptr<PrimaryObservable>> &obs, const std::vector<double> &expected) {
                    return obs.size() == 1 && obs[0]->getValues() == expected;
                },
                [](std::ostream &os, const std::vector<double> &expected) {
                    os << " == {";
                    std::copy(expected.begin(), expected.end(), std::ostream_iterator<double>(os, ", "));
                    os << "}";
                },
                values
        );
    }
}

TEST_CASE("EigenstateObservables: single set") {
    trompeloeil::sequence seq1;
    auto primary = std::make_shared<PrimaryObservableMock>();
    ALLOW_CALL(*primary, getHeader()).RETURN(std::vector<std::string>{"p"});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{1, 0, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{1});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{0, 1, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{2});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{0, 0, 1}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{3});

    trompeloeil::sequence seq2;
    auto secondary = std::make_shared<SecondaryObservableMock>();
    ALLOW_CALL(*secondary, getHeader()).RETURN(std::vector<std::string>{"s"});
    REQUIRE_CALL(*secondary, calculateForObservables(first_observable_values_eq({1}))).IN_SEQUENCE(seq2);
    REQUIRE_CALL(*secondary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq2).RETURN(std::vector<double>{4});
    REQUIRE_CALL(*secondary, calculateForObservables(first_observable_values_eq({2}))).IN_SEQUENCE(seq2);
    REQUIRE_CALL(*secondary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq2).RETURN(std::vector<double>{5});
    REQUIRE_CALL(*secondary, calculateForObservables(first_observable_values_eq({3}))).IN_SEQUENCE(seq2);
    REQUIRE_CALL(*secondary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq2).RETURN(std::vector<double>{6});

    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    EigenstateObservables eigenstateObservables(2, {primary}, {secondary}, {primary, secondary});
    eigenstateObservables.analyze(Eigensystem({0, 0.3, 1},
                                              {{1, 0, 0},
                                               {0, 1, 0},
                                               {0, 0, 1}}), logger);

    std::ostringstream out;
    eigenstateObservables.storeResult(out);
    REQUIRE(out.str() == "binStart p dp s ds \n"
                         "0 1.5 0 4.5 0 \n"
                         "0.5 3 0 6 0 \n");
}

TEST_CASE("EigenstateObservables: two sets - averaging") {
    trompeloeil::sequence seq1;
    auto primary = std::make_shared<PrimaryObservableMock>();
    ALLOW_CALL(*primary, getHeader()).RETURN(std::vector<std::string>{"p"});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{1, 0, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{1});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{0, 1, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{2});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{0, 0, 1}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{3});

    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{1, 0, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{4});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{0, 1, 0}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{5});
    REQUIRE_CALL(*primary, calculateForState(arma_eq(arma::cx_vec{0, 0, 1}))).IN_SEQUENCE(seq1);
    REQUIRE_CALL(*primary, getValues()).TIMES(AT_LEAST(1)).IN_SEQUENCE(seq1).RETURN(std::vector<double>{6});

    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    EigenstateObservables eigenstateObservables(2, {primary}, {}, {primary});
    eigenstateObservables.analyze(Eigensystem({0, 0.3, 1},
                                              {{1, 0, 0},
                                               {0, 1, 0},
                                               {0, 0, 1}}), logger);
    eigenstateObservables.analyze(Eigensystem({0, 0.7, 1},
                                              {{1, 0, 0},
                                               {0, 1, 0},
                                               {0, 0, 1}}), logger);

    std::ostringstream out;
    eigenstateObservables.storeResult(out);
    REQUIRE(out.str() == "binStart p dp \n"
                         "0 2.750 1.250 \n"
                         "0.5 4.250 1.250 \n");
}