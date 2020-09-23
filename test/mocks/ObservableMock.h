//
// Created by Piotr Kubala on 22/09/2020.
//

#ifndef MBL_ED_OBSERVABLEMOCK_H
#define MBL_ED_OBSERVABLEMOCK_H

#include <catch2/trompeloeil.hpp>

#include "evolution/Observable.h"

class ObservableMock : public trompeloeil::mock_interface<Observable> {
public:
    IMPLEMENT_CONST_MOCK0(getHeader);
    IMPLEMENT_CONST_MOCK0(getValues);
};

#endif //MBL_ED_OBSERVABLEMOCK_H
