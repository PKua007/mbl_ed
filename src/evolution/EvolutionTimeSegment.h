//
// Created by pkua on 29.05.2020.
//

#ifndef MBL_ED_EVOLUTIONTIMESEGMENT_H
#define MBL_ED_EVOLUTIONTIMESEGMENT_H

#include <iosfwd>

#include "utils/Assertions.h"

struct EvolutionTimeSegment {
    double maxTime{};
    std::size_t numSteps{};

    EvolutionTimeSegment() = default;
    EvolutionTimeSegment(double maxTime, std::size_t numSteps) : maxTime{maxTime}, numSteps{numSteps} {
        Expects(maxTime > 0);
        Expects(numSteps > 0);
    }

    friend std::istream &operator>>(std::istream &in, EvolutionTimeSegment &segment) {
        double maxTime_{};
        std::size_t numSteps_{};
        in >> maxTime_ >> numSteps_;
        if (!in)
            return in;

        segment = {maxTime_, numSteps_};
        return in;
    }
};


#endif //MBL_ED_EVOLUTIONTIMESEGMENT_H
