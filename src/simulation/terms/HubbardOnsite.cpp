//
// Created by Piotr Kubala on 09/02/2020.
//

#include <numeric>

#include "HubbardOnsite.h"

double HubbardOnsite::calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) {
    static_cast<void>(generator);

    auto bosonAccumulator = [](auto sum, auto numberOfParticles) {
        return sum + numberOfParticles*(numberOfParticles - 1);
    };
    return this->U / 2 * std::accumulate(vector.begin(), vector.end(), 0., bosonAccumulator);
}

HubbardOnsite::HubbardOnsite(double U) : U{U} {
    Expects(U >= 0);
}
