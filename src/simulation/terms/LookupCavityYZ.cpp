//
// Created by Piotr Kubala on 10/02/2020.
//

#include <numeric>

#include "LookupCavityYZ.h"
#include "simulation/HamiltonianGenerator.h"

void LookupCavityYZ::changeRealisation(std::size_t index) {
    Expects(index < this->cavityConstants.size());
    this->currentRealisation = this->cavityConstants[index];
}

double LookupCavityYZ::calculate(const HopData &hopData, const HamiltonianGenerator &generator)
{
    Expects(generator.getSiteDistance(hopData.fromSite, hopData.toSite) == 1);
    if (generator.usingPBC())
        throw std::runtime_error("LookupCavityYZ: PBC not supported");
    Assert(hopData.fromVector.size() <= this->cavityConstants.getNumberOfSites());

    std::size_t smallerSite = std::min(hopData.fromSite, hopData.toSite);
    double fromZTerm = this->calculateZTerm(hopData.fromVector);
    double toZTerm = this->calculateZTerm(hopData.toVector);
    return -this->U1 / hopData.fromVector.size() * (fromZTerm + toZTerm)
           * this->currentRealisation.siteEntries[smallerSite].y;
}

double LookupCavityYZ::calculateZTerm(const FockBase::Vector &vector) const {
    std::size_t siteIndex{};
    auto siteAccumulator = [&siteIndex, this](auto sum, auto element) {
        return sum + this->currentRealisation.siteEntries[siteIndex++].wannier * element;
    };
    return std::accumulate(vector.begin(), vector.end(), 0., siteAccumulator);
}
