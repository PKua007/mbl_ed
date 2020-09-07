//
// Created by Piotr Kubala on 09/02/2020.
//

#include <cmath>
#include <numeric>

#include "CavityLongInteraction.h"
#include "core/HamiltonianGenerator.h"

double CavityLongInteraction::calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) {
    Expects(!generator.usingPBC());

    std::size_t elementIndex{};
    auto plusMinusAccumulator = [&elementIndex, this](auto sum, auto element) {
        return sum + std::cos(2*M_PI*this->beta*(elementIndex++) + this->phi0 + this->phi0Bias) * element;
    };
    double populationImbalance = std::accumulate(vector.begin(), vector.end(), 0., plusMinusAccumulator);
    return -this->U1 / vector.size() * populationImbalance * populationImbalance;
}

CavityLongInteraction::CavityLongInteraction(double U1, double beta, double phi0, double phi0Bias)
        : U1{U1}, beta{beta}, phi0{phi0}, phi0Bias{phi0Bias}
{
    Expects(beta >= 0);
}

void CavityLongInteraction::setPhi0(double phi0_) {
    this->phi0 = phi0_;
}
