//
// Created by Piotr Kubala on 08/02/2020.
//

#include <catch2/catch.hpp>

#include "mocks/RNDMock.h"

#include "simulation/AveragingModels.h"

/*namespace {
    class MockHamiltonianGenerator {
    public:
        MAKE_MOCK1(resampleOnsiteEnergies, void(RND&));
        MAKE_MOCK1(setPhi0, void(double));
    };
}

TEST_CASE("OnsiteDisorderAveragingModel") {
    RNDMock rnd;
    MockHamiltonianGenerator generator;
    REQUIRE_CALL(generator, resampleOnsiteEnergies(trompeloeil::_))
        .LR_WITH(&_1 == &rnd);

    OnsiteDisorderAveragingModel<MockHamiltonianGenerator>::setupHamiltonianGenerator(generator, rnd, 0, 0);
}

TEST_CASE("UniformPhi0AveragingModel") {
    RNDMock rnd;
    MockHamiltonianGenerator generator;
    REQUIRE_CALL(generator, resampleOnsiteEnergies(trompeloeil::_))
         .LR_WITH(&_1 == &rnd);
    REQUIRE_CALL(generator, setPhi0(trompeloeil::_))
         .WITH(_1 == Approx(M_PI/3));

    UniformPhi0AveragingModel<MockHamiltonianGenerator>::setupHamiltonianGenerator(generator, rnd, 1, 3);
}

TEST_CASE("RandomPhi0AveragingModel") {
    RNDMock rnd;
    REQUIRE_CALL(rnd, getDouble())
        .RETURN(0.5);
    MockHamiltonianGenerator generator;
    REQUIRE_CALL(generator, resampleOnsiteEnergies(trompeloeil::_))
        .LR_WITH(&_1 == &rnd);
    REQUIRE_CALL(generator, setPhi0(trompeloeil::_))
        .WITH(_1 == Approx(M_PI));

    RandomPhi0AveragingModel<MockHamiltonianGenerator>::setupHamiltonianGenerator(generator, rnd, 0, 1);
}*/