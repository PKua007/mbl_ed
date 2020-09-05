//
// Created by Piotr Kubala on 16/03/2020.
//

#include <algorithm>
#include <numeric>
#include <cmath>

#include "utils/Assertions.h"
#include "QuasiperiodicDisorder.h"

QuasiperiodicDisorder::QuasiperiodicDisorder(double W, double beta, double phi0) : W{W}, beta{beta}, phi0{phi0} {
    Expects(W >= 0);
    Expects(beta >= 0);
}

double QuasiperiodicDisorder::calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) {
    static_cast<void>(generator);

    double energy{};
    for (std::size_t i{}; i < vector.size(); i++)
        energy += vector[i] * this->W * std::cos(2*M_PI*this->beta*i + this->phi0);

    return energy;
}

void QuasiperiodicDisorder::setPhi0(double phi0) {
    this->phi0 = phi0;
}