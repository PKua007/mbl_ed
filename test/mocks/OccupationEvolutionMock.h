//
// Created by Piotr Kubala on 22/09/2020.
//

#ifndef MBL_ED_OCCUPATIONEVOLUTIONMOCK_H
#define MBL_ED_OCCUPATIONEVOLUTIONMOCK_H

#include <catch2/trompeloeil.hpp>

#include "evolution/OservablesTimeEvolution.h"

class OccupationEvolutionMock : public OservablesTimeEvolution {
public:
    MAKE_MOCK4(perform, std::vector<TimeEvolutionEntry>(const std::vector<EvolutionTimeSegment> &,
                                                           const arma::cx_vec &, Evolver &, Logger &),
               override);
};

#endif //MBL_ED_OCCUPATIONEVOLUTIONMOCK_H
