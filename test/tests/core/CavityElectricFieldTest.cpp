//
// Created by Piotr Kubala on 07/12/2020.
//

#include <catch2/catch.hpp>

#include "core/observables/CavityElectricField.h"
#include "core/observables/CavityOnsiteOccupations.h"

namespace {
    class MockCavityOnsiteOccupations : public CavityOnsiteOccupations {
    public:
        [[nodiscard]] std::vector<double> getValues() const override {
            return {4, 5};
        }
    };
}

TEST_CASE("CavityElectricField") {
    CavityElectricField cavityElectricField;

    SECTION("header") {
        REQUIRE(cavityElectricField.getHeader() == std::vector<std::string>{"a"});
    }

    SECTION("values") {
        cavityElectricField.calculateForObservables({std::make_shared<MockCavityOnsiteOccupations>()});

        REQUIRE(cavityElectricField.getValues() == std::vector<double>{9});
    }
}