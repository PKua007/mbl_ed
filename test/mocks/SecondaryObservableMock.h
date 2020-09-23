//
// Created by Piotr Kubala on 22/09/2020.
//

#ifndef MBL_ED_SECONDARYOBSERVABLEMOCK_H
#define MBL_ED_SECONDARYOBSERVABLEMOCK_H

#include <catch2/trompeloeil.hpp>

#include "evolution/SecondaryObservable.h"

class SecondaryObservableMock : public trompeloeil::mock_interface<SecondaryObservable> {
public:
    IMPLEMENT_CONST_MOCK0(getHeader);
    IMPLEMENT_CONST_MOCK0(getValues);
    IMPLEMENT_MOCK1(calculateForObservables);
};

#endif //MBL_ED_SECONDARYOBSERVABLEMOCK_H
