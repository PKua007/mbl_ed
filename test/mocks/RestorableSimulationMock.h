//
// Created by Piotr Kubala on 01/09/2020.
//

#ifndef MBL_ED_RESTORABLESIMULATIONMOCK_H
#define MBL_ED_RESTORABLESIMULATIONMOCK_H

#include <catch2/trompeloeil.hpp>

#include "simulation/RestorableSimulation.h"

class RestorableSimulationMock : public trompeloeil::mock_interface<RestorableSimulation> {
    IMPLEMENT_CONST_MOCK1(storeState);
    IMPLEMENT_MOCK1(joinRestoredState);
    IMPLEMENT_MOCK0(clear);
    IMPLEMENT_MOCK1(seedRandomGenerators);
    IMPLEMENT_MOCK3(performSimulation);
    IMPLEMENT_CONST_MOCK0(getTagName);
};

#endif //MBL_ED_RESTORABLESIMULATIONMOCK_H
