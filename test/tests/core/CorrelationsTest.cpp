//
// Created by Piotr Kubala on 23/09/2020.
//

#include <catch2/catch.hpp>

#include "core/observables/Correlations.h"
#include "core/observables/OnsiteOccupations.h"
#include "core/observables/OnsiteOccupationsSquared.h"

namespace {
    class MockOnsiteOccupations : public OnsiteOccupations {
    public:
        std::vector<double> getValues() const override {
            return {4, 5, 6, 7, 8};
        }
    };

    class MockOnsiteOccupationsSquared : public OnsiteOccupationsSquared {
    private:
        SymmetricMatrix<double> r;
    public:
        MockOnsiteOccupationsSquared() : r(5) {
            r(0, 0) = 1;
            r(1, 0) = 2; r(1, 1) = 3;
            r(2, 0) = 4; r(2, 1) = 5; r(2, 2) = 6;
            r(3, 0) = 7; r(3, 1) = 8; r(3, 2) = 9; r(3, 3) = 10;
            r(4, 0) = 11; r(4, 1) = 12; r(4, 2) = 13; r(4, 3) = 14; r(4, 4) = 14;
        }

        const SymmetricMatrix<double> &getOccupationsSquared() const override {
            return this->r;
        }
    };
}

TEST_CASE("Correlations") {
    Correlations correlations(5, 1);

    SECTION("header") {
        REQUIRE(correlations.getHeader() == std::vector<std::string>{"Gm1_1", "Gm1_2"});
    }

    SECTION("values") {
        correlations.calculateForObservables({std::make_shared<MockOnsiteOccupations>(),
                                              std::make_shared<MockOnsiteOccupationsSquared>()});

        REQUIRE(correlations.getValues() == std::vector<double>{-29, -27});
    }
}