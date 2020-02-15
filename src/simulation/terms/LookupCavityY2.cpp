//
// Created by Piotr Kubala on 14/02/2020.
//

#include "LookupCavityY2.h"
#include "simulation/HamiltonianGenerator.h"

double LookupCavityY2::calculate(const HopData &firstHopData, const HopData &secondHopData,
                                 const HamiltonianGenerator &generator) {
    Expects(generator.getSiteDistance(firstHopData.fromSite, firstHopData.toSite) == 1);
    Expects(generator.getSiteDistance(secondHopData.fromSite, secondHopData.toSite) == 1);
    if (generator.usingPBC())
        throw std::runtime_error("LookupCavityY2: PBC not supported");
    std::size_t numberOfSites = firstHopData.fromVector.size();
    Assert(numberOfSites <= this->cavityConstants.getNumberOfSites());

    std::size_t firstSmallerSite = std::min(firstHopData.fromSite, firstHopData.toSite);
    std::size_t secondSmallerSite = std::min(secondHopData.fromSite, secondHopData.toSite);
    return -this->U1 / numberOfSites * this->currentRealisation.siteEntries[firstSmallerSite].y
           * this->currentRealisation.siteEntries[secondSmallerSite].y;
}

void LookupCavityY2::changeRealisation(std::size_t index) {
    Expects(index < this->cavityConstants.size());
    this->currentRealisation = this->cavityConstants[index];
}