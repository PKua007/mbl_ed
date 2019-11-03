//
// Created by pkua on 01.11.2019.
//

#include <catch2/catch.hpp>

#include "FockBaseGenerator.h"
#include "HamiltonianGenerator.h"
#include "Assertions.h"

namespace {
    class MockPeriodicHamiltonianGenerator : public HamiltonianGenerator {
    public:
        mutable bool periodicUsed = false;

        explicit MockPeriodicHamiltonianGenerator(const FockBase &fockBase) : HamiltonianGenerator(fockBase) { }

        [[nodiscard]] double getDiagonalElement(const FockBase::Vector &vector) const override {
            return *(this->fockBase.findIndex(vector));
        }

        [[nodiscard]] double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const override {
            Expects(fromSiteIndex < this->fockBase.getNumberOfSites());
            Expects(toSiteIndex < this->fockBase.getNumberOfSites());
            Expects(std::abs(static_cast<int>(toSiteIndex - fromSiteIndex)) == 1
                    || std::abs(static_cast<int>(toSiteIndex - fromSiteIndex))
                       == static_cast<int>(this->fockBase.getNumberOfSites() - 1));
            if (std::abs(static_cast<int>(toSiteIndex) - static_cast<int>(fromSiteIndex)) == 1) {
                return -1;
            } else {
                this->periodicUsed = true;
                return 0;
            }
        }
    };
}

namespace {
    class MockNonPeriodicHamiltonianGenerator : public HamiltonianGenerator {
    public:
        explicit MockNonPeriodicHamiltonianGenerator(const FockBase &fockBase) : HamiltonianGenerator(fockBase, false) {

        }

        [[nodiscard]] double getDiagonalElement(const FockBase::Vector &vector) const override {
            return *(this->fockBase.findIndex(vector));
        }

        [[nodiscard]] double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const override {
            Expects(fromSiteIndex < this->fockBase.getNumberOfSites());
            Expects(toSiteIndex < this->fockBase.getNumberOfSites());
            Expects(std::abs(static_cast<int>(toSiteIndex - fromSiteIndex)) == 1);

            return -1;
        }
    };
}

TEST_CASE("HamiltonianGenerator: 2 bosons in 3 sites") {
    FockBaseGenerator baseGenerator;
    FockBase fockBase = baseGenerator.generate(3, 2);
    MockPeriodicHamiltonianGenerator hamiltonianGenerator(fockBase);

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

TEST_CASE("PBC") {
    SECTION("Periodic BC - periodic hopping should be present") {
        FockBaseGenerator baseGenerator;
        FockBase fockBase = baseGenerator.generate(3, 2);
        MockPeriodicHamiltonianGenerator hamiltonianGenerator(fockBase);

        static_cast<void>(hamiltonianGenerator.generate());

        REQUIRE(hamiltonianGenerator.periodicUsed);
    }

    SECTION("Non periodic BC - periodic hopping should not be present") {
        FockBaseGenerator baseGenerator;
        FockBase fockBase = baseGenerator.generate(3, 2);
        MockNonPeriodicHamiltonianGenerator hamiltonianGenerator(fockBase);

        REQUIRE_NOTHROW(hamiltonianGenerator.generate());
    }
}