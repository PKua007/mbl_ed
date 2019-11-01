//
// Created by pkua on 01.11.2019.
//

#include <catch2/catch.hpp>

#include "FockBaseGenerator.h"
#include "HamiltonianGenerator.h"
#include "Assertions.h"

namespace {
    class MockHamiltonianGenerator : public HamiltonianGenerator {
    public:
        explicit MockHamiltonianGenerator(const FockBase &fockBase) : HamiltonianGenerator(fockBase) { }

        [[nodiscard]] double getDiagonalElement(const FockBase::Vector &vector) const override {
            return *(this->fockBase.findIndex(vector));
        }

        [[nodiscard]] double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const override {
            Expects(fromSiteIndex < this->fockBase.getNumberOfSites());
            Expects(toSiteIndex < this->fockBase.getNumberOfSites());
            Expects(toSiteIndex - fromSiteIndex == 1);
            return -1;
        }
    };
}

TEST_CASE("HamiltonianGenerator: 2 bosons in 3 sites") {
    FockBaseGenerator baseGenerator;
    FockBase fockBase = baseGenerator.generate(3, 2);
    MockHamiltonianGenerator hamiltonianGenerator(fockBase);

    arma::mat result = hamiltonianGenerator.generate();

    arma::mat expected = {{ 0,       -M_SQRT2,  0,  0,        0,        0},
                          {-M_SQRT2,  1,       -1, -M_SQRT2,  0,        0},
                          { 0,       -1,        2,  0,       -1,        0},
                          { 0,       -M_SQRT2,  0,  3,       -M_SQRT2,  0},
                          { 0,        0,       -1, -M_SQRT2,  4,       -M_SQRT2},
                          { 0,        0,        0,  0,       -M_SQRT2,  5}};
    REQUIRE(arma::any(arma::vectorise(result - expected) > -0.000001));
    REQUIRE(arma::any(arma::vectorise(result - expected) < 0.000001));
}