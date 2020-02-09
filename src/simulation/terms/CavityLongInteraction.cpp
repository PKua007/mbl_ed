//
// Created by Piotr Kubala on 09/02/2020.
//

#include <cmath>
#include <numeric>

#include "CavityLongInteraction.h"

double CavityLongInteraction::calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) {
    static_cast<void>(generator);

    std::size_t elementIndex{};
    auto plusMinusAccumulator = [&elementIndex, this](auto sum, auto element) {
        return sum + std::cos(beta*(elementIndex++) + phi0) * element;
    };
    double populationImbalance = std::accumulate(vector.begin(), vector.end(), 0., plusMinusAccumulator);
    return -this->U1 / vector.size() * populationImbalance * populationImbalance;
}

CavityLongInteraction::CavityLongInteraction(double U1, double beta, double phi0) : U1{U1}, beta{beta}, phi0{phi0} {
    Expects(U1 >= 0);
    Expects(beta >= 0);
}

void CavityLongInteraction::setPhi0(double phi0) {
    this->phi0 = phi0;
}
