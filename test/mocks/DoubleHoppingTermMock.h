//
// Created by Piotr Kubala on 14/02/2020.
//

#ifndef MBL_ED_DOUBLEHOPPINGTERMMOCK_H
#define MBL_ED_DOUBLEHOPPINGTERMMOCK_H

#include <catch2/trompeloeil.hpp>

#include "simulation/DoubleHoppingTerm.h"

class DoubleHoppingTermMock : public trompeloeil::mock_interface<DoubleHoppingTerm> {
    IMPLEMENT_MOCK3(calculate);
};

#endif //MBL_ED_DOUBLEHOPPINGTERMMOCK_H
