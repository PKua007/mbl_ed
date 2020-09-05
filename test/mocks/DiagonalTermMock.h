//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_DIAGONALTERMMOCK_H
#define MBL_ED_DIAGONALTERMMOCK_H

#include <catch2/trompeloeil.hpp>

#include "core/DiagonalTerm.h"

class DiagonalTermMock : public trompeloeil::mock_interface<DiagonalTerm> {
    IMPLEMENT_MOCK2(calculate);
};

#endif //MBL_ED_DIAGONALTERMMOCK_H
