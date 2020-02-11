//
// Created by Piotr Kubala on 10/02/2020.
//

#include <numeric>

#include "LookupCavityZ2.h"
#include "simulation/HamiltonianGenerator.h"

double LookupCavityZ2::calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) {
    Expects(!generator.usingPBC());
    Expects(vector.size() <= this->cavityConstants.getNumberOfSites());

    std::size_t siteIndex{};
    auto siteAccumulator = [&siteIndex, this](auto sum, auto element) {
        return sum + this->currentRealisation.siteEntries[siteIndex++].wannier * element;
    };
    double populationImbalance = std::accumulate(vector.begin(), vector.end(), 0., siteAccumulator);
    return -this->U1 / vector.size() * populationImbalance * populationImbalance;
}

void LookupCavityZ2::changeRealisation(std::size_t index) {
    Expects(index < this->cavityConstants.size());
    this->currentRealisation = this->cavityConstants[index];
}
