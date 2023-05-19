//
// Created by Piotr Kubala on 03/10/2020.
//

#ifndef MBL_ED_FIRSTOBSERVABLEVALUESTROMPELOEILMATCHER_H
#define MBL_ED_FIRSTOBSERVABLEVALUESTROMPELOEILMATCHER_H

#include <iterator>

#include <catch2/trompeloeil.hpp>

#include "core/PrimaryObservable.h"

inline auto first_observable_values_eq(const std::vector<double> &values) {
    return trompeloeil::make_matcher<std::vector<std::shared_ptr<PrimaryObservable>>>(
            [](const std::vector<std::shared_ptr<PrimaryObservable>> &obs, const std::vector<double> &expected) {
                return obs.size() == 1 && obs[0]->getValues() == expected;
            },
            [](std::ostream &os, const std::vector<double> &expected) {
                os << " == {";
                std::copy(expected.begin(), expected.end(), std::ostream_iterator<double>(os, ", "));
                os << "}";
            },
            values
    );
}

#endif //MBL_ED_FIRSTOBSERVABLEVALUESTROMPELOEILMATCHER_H
