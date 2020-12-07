//
// Created by Piotr Kubala on 23/09/2020.
//

#include <catch2/catch.hpp>

#include "core/observables/OnsiteFluctuations.h"
#include "core/observables/OnsiteOccupations.h"
#include "core/observables/OnsiteOccupationsSquared.h"

namespace {
    class MockOnsiteOccupations : public OnsiteOccupations {
    public:
        std::vector<double> getValues() const override {
            return {4, 5};
        }
    };

    class MockOnsiteOccupationsSquared : public OnsiteOccupationsSquared {
    private:
        SymmetricMatrix<double> r;
    public:
        MockOnsiteOccupationsSquared() : r(2) {
            r(0, 0) = 24;
            r(1, 0) = 25; r(1, 1) = 26;
        }

        const SymmetricMatrix<double> &getOccupationsSquared() const override {
            return this->r;
        }
    };
}

TEST_CASE("OnsiteFluctuations") {
    OnsiteFluctuations onsiteFluctuations(2);

    SECTION("header") {
        REQUIRE(onsiteFluctuations.getHeader() == std::vector<std::string>{"rho_1", "rho_2"});
    }

    SECTION("values") {
        onsiteFluctuations.calculateForObservables({std::make_shared<MockOnsiteOccupations>(),
                                                    std::make_shared<MockOnsiteOccupationsSquared>()});

        REQUIRE(onsiteFluctuations.getValues() == std::vector<double>{8, 1});
    }
}