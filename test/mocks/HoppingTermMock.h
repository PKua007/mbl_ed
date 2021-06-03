//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_HOPPINGTERMMOCK_H
#define MBL_ED_HOPPINGTERMMOCK_H

#include <catch2/trompeloeil.hpp>

#include "core/HoppingTerm.h"

class HoppingTermMock : public trompeloeil::mock_interface<HoppingTerm> {
    IMPLEMENT_CONST_MOCK2(calculate);
    IMPLEMENT_CONST_MOCK0(getHoppingDistances);
};

#endif //MBL_ED_HOPPINGTERMMOCK_H
