//
// Created by Piotr Kubala on 21/09/2020.
//

#include "Correlations.h"
#include "utils/Assertions.h"
#include "evolution/SymmetricMatrix.h"
#include "OnsiteOccupationsSquared.h"
#include "OnsiteOccupations.h"

std::vector<std::string> Correlations::getHeader() const {
    std::vector<std::string> headerStrings;
    for (std::size_t d = 1; d < this->borderlessNumOfSites; d++)
        headerStrings.push_back("Gm" + std::to_string(this->marginSize) + "_" + std::to_string(d));
    return headerStrings;
}

std::vector<double> Correlations::getValues() const {
    return this->G_d;
}

Correlations::Correlations(std::size_t numOfSites, std::size_t marginSize)
        : numOfSites{numOfSites}, marginSize{marginSize}, borderlessNumOfSites{numOfSites - 2*marginSize}
{
    Expects(numOfSites > 0);
    Expects(numOfSites - 2*marginSize >= 2);
    this->G_d.resize(this->borderlessNumOfSites - 1);
}

void Correlations::calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) {
    std::fill(this->G_d.begin(), this->G_d.end(), 0.);

    auto numParticles = findObservable<OnsiteOccupations>(primaryObservables).getValues();
    auto numParticlesSquared = findObservable<OnsiteOccupationsSquared>(primaryObservables).getOccupationsSquared();
    Assert(numParticles.size() == this->numOfSites);
    Assert(numParticlesSquared.size() == this->numOfSites);

    for (std::size_t d = 1; d < this->borderlessNumOfSites; d++) {
        std::size_t minFirstSite = this->marginSize;
        std::size_t maxFirstSite = this->numOfSites - d - this->marginSize - 1;
        for (std::size_t firstSite = minFirstSite; firstSite <= maxFirstSite; firstSite++) {
            std::size_t secondSite = firstSite + d;
            double singleG = numParticlesSquared(firstSite, secondSite) -
                             numParticles[firstSite] * numParticles[secondSite];
            this->G_d[d - 1] += singleG;
        }
        std::size_t numSummands = this->numOfSites - 2*this->marginSize - d;
        this->G_d[d - 1] /= numSummands;
    }
}