//
// Created by Piotr Kubala on 22/09/2020.
//

#ifndef MBL_ED_EVOLVERMOCK_H
#define MBL_ED_EVOLVERMOCK_H

#include <catch2/trompeloeil.hpp>

#include "evolution/Evolver.h"

class EvolverMock : public trompeloeil::mock_interface<Evolver> {
public:
    IMPLEMENT_MOCK3(prepareFor);
    IMPLEMENT_MOCK0(evolve);
    IMPLEMENT_CONST_MOCK0(getCurrentState);
    IMPLEMENT_CONST_MOCK0(getDt);
};

#endif //MBL_ED_EVOLVERMOCK_H