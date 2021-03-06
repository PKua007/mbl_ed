//
// Created by Piotr Kubala on 21/09/2020.
//

#include "OnsiteFluctuations.h"
#include "utils/Assertions.h"
#include "evolution/SymmetricMatrix.h"
#include "OnsiteOccupationsSquared.h"
#include "OnsiteOccupations.h"

std::vector<std::string> OnsiteFluctuations::getHeader() const {
    std::vector<std::string> headerStrings;
    for (std::size_t i = 0; i < numOfSites; i++)
        headerStrings.push_back("rho_" + std::to_string(i + 1));
    return headerStrings;
}

std::vector<double> OnsiteFluctuations::getValues() const {
    return this->rho_i;
}

OnsiteFluctuations::OnsiteFluctuations(std::size_t numOfSites) : numOfSites{numOfSites}, rho_i(numOfSites) {
    Expects(numOfSites > 0);
}

void
OnsiteFluctuations::calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables)
{
    auto numParticles = findObservable<OnsiteOccupations>(primaryObservables).getValues();
    auto numParticlesSquared = findObservable<OnsiteOccupationsSquared>(primaryObservables).getOccupationsSquared();
    Assert(numParticles.size() == this->numOfSites);
    Assert(numParticlesSquared.size() == this->numOfSites);

    for (std::size_t i = 0; i < this->numOfSites; i++)
        this->rho_i[i] = numParticlesSquared(i, i) - std::pow(numParticles[i], 2);
}