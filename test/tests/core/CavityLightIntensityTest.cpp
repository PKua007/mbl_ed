//
// Created by Piotr Kubala on 07/12/2020.
//

#include <catch2/catch.hpp>

#include "core/observables/CavityLightIntensity.h"
#include "core/observables/CavityOnsiteOccupationsSquared.h"

namespace {
    class MockCavityOnsiteOccupationsSquared : public CavityOnsiteOccupationsSquared {
    private:
        SymmetricMatrix<double> r;
    public:
        MockCavityOnsiteOccupationsSquared() : r(2) {
            r(0, 0) = 1;
            r(1, 0) = 2; r(1, 1) = 3;
        }

        const SymmetricMatrix<double> &getOccupationsSquared() const override {
            return this->r;
        }
    };
}

TEST_CASE("CavityLightIntensity") {
    CavityLightIntensity lightIntensity;

    SECTION("header") {
        REQUIRE(lightIntensity.getHeader() == std::vector<std::string>{"ada"});
    }

    SECTION("values") {
        lightIntensity.calculateForObservables({std::make_shared<MockCavityOnsiteOccupationsSquared>()});

        REQUIRE(lightIntensity.getValues() == std::vector<double>{8});
    }
}