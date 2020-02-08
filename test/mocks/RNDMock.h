//
// Created by Piotr Kubala on 08/02/2020.
//

#ifndef MBL_ED_RNDMOCK_H
#define MBL_ED_RNDMOCK_H

#include <catch2/trompeloeil.hpp>

#include "simulation/RND.h"

class RNDMock : public RND {
public:
    MAKE_MOCK0(getDouble, double(), override);
};

#endif //MBL_ED_RNDMOCK_H
