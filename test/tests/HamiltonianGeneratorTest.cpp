//
// Created by pkua on 01.11.2019.
//

#include <catch2/catch.hpp>

#include "matchers/ArmaApproxEqualMatcher.h"
#include "simulation/FockBaseGenerator.h"
#include "simulation/HamiltonianGenerator.h"
#include "utils/Assertions.h"

namespace {
    class MockPeriodicHamiltonianGenerator : public HamiltonianGenerator {
    public:
        mutable bool periodicUsed = false;

        explicit MockPeriodicHamiltonianGenerator(std::unique_ptr<FockBase> fockBase)
                : HamiltonianGenerator(std::move(fockBase), true)
        { }

        [[nodiscard]] double getDiagonalElement(const FockBase::Vector &vector) const override {
            return *(this->fockBase->findIndex(vector));
        }

        [[nodiscard]] double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const override {
            Expects(fromSiteIndex < this->fockBase->getNumberOfSites());
            Expects(toSiteIndex < this->fockBase->getNumberOfSites());
            Expects(std::abs(static_cast<int>(toSiteIndex - fromSiteIndex)) == 1
                    || std::abs(static_cast<int>(toSiteIndex - fromSiteIndex))
                       == static_cast<int>(this->fockBase->getNumberOfSites() - 1));
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
        explicit MockNonPeriodicHamiltonianGenerator(std::unique_ptr<FockBase> fockBase)
                : HamiltonianGenerator(std::move(fockBase), false)
        { }

        [[nodiscard]] double getDiagonalElement(const FockBase::Vector &vector) const override {
            return *(this->fockBase->findIndex(vector));
        }

        [[nodiscard]] double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const override {
            Expects(fromSiteIndex < this->fockBase->getNumberOfSites());
            Expects(toSiteIndex < this->fockBase->getNumberOfSites());
            Expects(std::abs(static_cast<int>(toSiteIndex - fromSiteIndex)) == 1);

            return -1;
        }
    };
}

TEST_CASE("HamiltonianGenerator: 2 bosons in 3 sites") {
    FockBaseGenerator baseGenerator;
    auto fockBase = baseGenerator.generate(3, 2);
    MockPeriodicHamiltonianGenerator hamiltonianGenerator(std::move(fockBase));

    arma::mat result = hamiltonianGenerator.generate();

    arma::mat expected = {{ 0,       -M_SQRT2,  0,  0,        0,        0},
                          {-M_SQRT2,  1,       -1, -M_SQRT2,  0,        0},
                          { 0,       -1,        2,  0,       -1,        0},
                          { 0,       -M_SQRT2,  0,  3,       -M_SQRT2,  0},
                          { 0,        0,       -1, -M_SQRT2,  4,       -M_SQRT2},
                          { 0,        0,        0,  0,       -M_SQRT2,  5}};
    REQUIRE_THAT(result, IsApproxEqual(expected, 1e-8));
}

TEST_CASE("HamiltonianGenerator: PBC") {
    SECTION("Periodic BC - periodic hopping should be present") {
        FockBaseGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(3, 2);
        MockPeriodicHamiltonianGenerator hamiltonianGenerator(std::move(fockBase));

        static_cast<void>(hamiltonianGenerator.generate());

        REQUIRE(hamiltonianGenerator.periodicUsed);
    }

    SECTION("Non periodic BC - periodic hopping should not be present") {
        FockBaseGenerator baseGenerator;
        auto fockBase = baseGenerator.generate(3, 2);
        MockNonPeriodicHamiltonianGenerator hamiltonianGenerator(std::move(fockBase));

        REQUIRE_NOTHROW(hamiltonianGenerator.generate());
    }
}

TEST_CASE("HamiltonianGenerator: site distance") {
    FockBaseGenerator baseGenerator;
    auto evenBase = baseGenerator.generate(6, 1);
    auto oddBase = baseGenerator.generate(7, 1);

    SECTION("Periodic BC - site distance calculated normally") {
        SECTION("even size") {
            MockPeriodicHamiltonianGenerator hamiltonianGenerator(std::move(evenBase));

            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 1) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(4, 5) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 3) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 4) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 5) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 5) == 1);
        }

        SECTION("odd size") {
            MockPeriodicHamiltonianGenerator hamiltonianGenerator(std::move(oddBase));

            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 1) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(5, 6) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 3) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 4) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 5) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 6) == 1);
        }
    }

    SECTION("Non-periodic BC - site distance calculated periodically") {
        SECTION("even size") {
            MockNonPeriodicHamiltonianGenerator hamiltonianGenerator(std::move(evenBase));

            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 1) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(4, 5) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 3) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 4) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 5) == 4);
            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 5) == 5);
        }

        SECTION("odd size") {
            MockNonPeriodicHamiltonianGenerator hamiltonianGenerator(std::move(oddBase));

            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 1) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(5, 6) == 1);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 3) == 2);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 4) == 3);
            REQUIRE(hamiltonianGenerator.getSiteDistance(1, 5) == 4);
            REQUIRE(hamiltonianGenerator.getSiteDistance(0, 6) == 6);
        }
    }
}