//
// Created by Piotr Kubala on 20/03/2020.
//

#include "ListOnsite.h"

#include "utils/Assertions.h"

double ListOnsite::calculate(const FockBasis::Vector &vector, const HamiltonianGenerator &generator) {
    Expects(vector.size() == this->onsitePotential.size());

    static_cast<void>(generator);

    double energy{};
    for (std::size_t i{}; i < vector.size(); i++)
        energy += vector[i] * this->onsitePotential[i];

    return energy;
}
