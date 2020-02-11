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

double LookupCavityYZ::calculate(const FockBase::Vector &from, const FockBase::Vector &to, std::size_t fromSite,
                                 std::size_t toSite, const HamiltonianGenerator &generator)
{
    Expects(generator.getSiteDistance(fromSite, toSite) == 1);
    if (generator.usingPBC())
        throw std::runtime_error("LookupCavityYZ: PBC not supported");
    Assert(from.size() <= this->cavityConstants.getNumberOfSites());

    std::size_t smallerSite = std::min(fromSite, toSite);
    double fromZTerm = this->calculateZTerm(from);
    double toZTerm = this->calculateZTerm(to);
    return -this->U1 / from.size() * (fromZTerm + toZTerm) * this->currentRealisation.siteEntries[smallerSite].y;
}

double LookupCavityYZ::calculateZTerm(const FockBase::Vector &vector) const {
    std::size_t siteIndex{};
    auto siteAccumulator = [&siteIndex, this](auto sum, auto element) {
        return sum + this->currentRealisation.siteEntries[siteIndex++].wannier * element;
    };
    return std::accumulate(vector.begin(), vector.end(), 0., siteAccumulator);
}
