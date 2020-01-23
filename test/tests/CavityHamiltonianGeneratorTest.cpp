//
// Created by pkua on 08.11.2019.
//

#include <catch2/catch.hpp>
#include <utility>
#include <simulation/FockBaseGenerator.h>

#include "matchers/ArmaApproxEqualMatcher.h"
#include "simulation/CavityHamiltonianGenerator.h"

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
    auto fockBase = baseGenerator.generate(3, 3);
    auto disorderGenerator = std::make_unique<SequenceDisorderGenerator>(std::vector<double>{1, 2, 3});
    CavityHamiltonianGeneratorParameters params{};
    CavityHamiltonianGenerator hamiltonianGenerator(std::move(fockBase), params, std::move(disorderGenerator));

    arma::mat result = hamiltonianGenerator.generate();

    arma::mat expected = arma::diagmat(arma::vec{3, 4, 5, 5, 6, 7, 6, 7, 8, 9});
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("CavityHamiltonianGenerator: short interactions") {
    FockBaseGenerator baseGenerator;
    auto fockBase = baseGenerator.generate(2, 4);
    CavityHamiltonianGeneratorParameters params{};
    params.U = 2;
    CavityHamiltonianGenerator hamiltonianGenerator(std::move(fockBase), params,
                                                    std::make_unique<ConstantDisorderGenerator>(0));

    arma::mat result = hamiltonianGenerator.generate();

    arma::mat expected = arma::diagmat(arma::vec{12, 6, 4, 6, 12});
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("CavityHamiltonianGenerator: long interactions") {
    SECTION("+-1 interactions") {
        FockBaseGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(3, 3);
        CavityHamiltonianGeneratorParameters params{};
        params.U1 = 3;
        params.beta = 0.5;
        CavityHamiltonianGenerator hamiltonianGenerator(std::move(fockBase), params,
                                                        std::make_unique<ConstantDisorderGenerator>(0));

        arma::mat result = hamiltonianGenerator.generate();

        arma::mat expected = arma::diagmat(arma::vec{-9, -1, -9, -1, -1, -9, -9, -1, -1, -9});
        REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
    }

    SECTION("cos interactions: beta=1/12, phi0=pi/6") {
        FockBaseGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(3, 2);
        CavityHamiltonianGeneratorParameters params{};
        params.U1 = 3;
        params.beta = 1./12;
        params.phi0 = M_PI/6;
        CavityHamiltonianGenerator hamiltonianGenerator(std::move(fockBase), params,
                                                        std::make_unique<ConstantDisorderGenerator>(0));

        arma::mat result = hamiltonianGenerator.generate();

        arma::mat expected = arma::diagmat(arma::vec{-3, -1-std::sqrt(3)/2, -3./4, -1, -1./4, 0});
        REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
    }
}