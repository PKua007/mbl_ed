//
// Created by pkua on 29.05.2020.
//

#ifndef MBL_ED_EVOLUTIONTIMESEGMENT_H
#define MBL_ED_EVOLUTIONTIMESEGMENT_H

#include <iosfwd>

#include "utils/Assertions.h"

/**
 * @brief A class representing time segment in time evolution.
 * @details For example, @a maxTime = 3, @a numSteps = 3 represents time interval od length 3 divided into 3 parts.
 */
struct EvolutionTimeSegment {
    double maxTime{};
    std::size_t numSteps{};

    EvolutionTimeSegment() = default;
    EvolutionTimeSegment(double maxTime, std::size_t numSteps) : maxTime{maxTime}, numSteps{numSteps} {
        Expects(maxTime > 0);
        Expects(numSteps > 0);
    }

    /**
     * @brief Stream extraction operator expecting format "maxTime numSteps"
     */
    friend std::istream &operator>>(std::istream &in, EvolutionTimeSegment &segment) {
        double maxTime_{};
        std::size_t numSteps_{};
        in >> maxTime_ >> numSteps_;
        if (!in)
            return in;

        segment = {maxTime_, numSteps_};
        return in;
    }

    friend bool operator==(const EvolutionTimeSegment &lhs, const EvolutionTimeSegment &rhs) {
        return std::tie(lhs.maxTime, lhs.numSteps) == std::tie(rhs.maxTime, rhs.numSteps);
    }

    friend bool operator!=(const EvolutionTimeSegment &lhs, const EvolutionTimeSegment &rhs) {
        return !(rhs == lhs);
    }
};


#endif //MBL_ED_EVOLUTIONTIMESEGMENT_H
