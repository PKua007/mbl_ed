//
// Created by Piotr Kubala on 22/09/2020.
//

#ifndef MBL_ED_PRIMARYOBSERVABLEMOCK_H
#define MBL_ED_PRIMARYOBSERVABLEMOCK_H

#include <catch2/trompeloeil.hpp>

#include "evolution/PrimaryObservable.h"

class PrimaryObservableMock : public trompeloeil::mock_interface<PrimaryObservable> {
public:
    IMPLEMENT_CONST_MOCK0(getHeader);
    IMPLEMENT_CONST_MOCK0(getValues);
    IMPLEMENT_MOCK1(calculateForState);
};

#endif //MBL_ED_PRIMARYOBSERVABLEMOCK_H
