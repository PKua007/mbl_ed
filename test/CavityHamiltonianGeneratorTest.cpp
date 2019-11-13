//
// Created by pkua on 08.11.2019.
//

#include <catch2/catch.hpp>
#include <utility>
#include <FockBaseGenerator.h>

#include "CavityHamiltonianGenerator.h"

namespace {
    class ConstantDisorderGenerator {
    private:
        double constant{};

    public:
        explicit ConstantDisorderGenerator(double constant) : constant(constant) { }

        double operator()() { return constant; }
    };

    class SequenceDisorderGenerator {
    private:
        std::size_t index{};
        std::vector<double> sequence;

    public:
        explicit SequenceDisorderGenerator(std::vector<double> sequence) : sequence(std::move(sequence)) { }

        double operator()() { return this->sequence.at(index++); }
    };
}

TEST_CASE("CavityHamiltonianGenerator: onsite energy") {
    FockBaseGenerator baseGenerator;
    FockBase fockBase = baseGenerator.generate(3, 3);
    auto disorderGenerator = std::make_unique<SequenceDisorderGenerator>(std::vector<double>{1, 2, 3});
    CavityHamiltonianGenerator hamiltonianGenerator(fockBase, 0, 0, 0, std::move(disorderGenerator));

    arma::mat result = hamiltonianGenerator.generate();

    arma::mat expected = arma::diagmat(arma::vec{3, 4, 5, 5, 6, 7, 6, 7, 8, 9});
    REQUIRE(arma::any(arma::vectorise(result - expected) > -0.000001));
    REQUIRE(arma::any(arma::vectorise(result - expected) < 0.000001));
}

TEST_CASE("CavityHamiltonianGenerator: short interactions") {
    FockBaseGenerator baseGenerator;
    FockBase fockBase = baseGenerator.generate(2, 4);
    CavityHamiltonianGenerator hamiltonianGenerator(fockBase, 0, 0, 1,
                                                    std::make_unique<ConstantDisorderGenerator>(0));

    arma::mat result = hamiltonianGenerator.generate();

    arma::mat expected = arma::diagmat(arma::vec{12, 6, 4, 6, 12});
    REQUIRE(arma::any(arma::vectorise(result - expected) > -0.000001));
    REQUIRE(arma::any(arma::vectorise(result - expected) < 0.000001));
}

TEST_CASE("CavityHamiltonianGenerator: long interactions") {
    FockBaseGenerator baseGenerator;
    FockBase fockBase = baseGenerator.generate(3, 3);
    CavityHamiltonianGenerator hamiltonianGenerator(fockBase, 0, 0, 0,
                                                    std::make_unique<ConstantDisorderGenerator>(0));

    arma::mat result = hamiltonianGenerator.generate();

    arma::mat expected = arma::diagmat(arma::vec{-9, -1, -9, -1, -1, -9, -9, -1, -1, -9});
    REQUIRE(arma::any(arma::vectorise(result - expected) > -0.000001));
    REQUIRE(arma::any(arma::vectorise(result - expected) < 0.000001));
}