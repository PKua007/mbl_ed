//
// Created by pkua on 29.05.2020.
//

#ifndef MBL_ED_EVOLUTIONTIMESEGMENT_H
#define MBL_ED_EVOLUTIONTIMESEGMENT_H

#include <cstddef>

#include "utils/Assertions.h"

struct EvolutionTimeSegment {
    const double maxTime{};
    const std::size_t numSteps{};

    EvolutionTimeSegment() = default;
    EvolutionTimeSegment(double maxTime, std::size_t numSteps) : maxTime{maxTime}, numSteps{numSteps} {
        Expects(maxTime > 0);
        Expects(numSteps > 0);
    }
};


#endif //MBL_ED_EVOLUTIONTIMESEGMENT_H
