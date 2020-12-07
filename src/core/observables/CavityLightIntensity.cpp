//
// Created by Piotr Kubala on 06/12/2020.
//

#include "CavityLightIntensity.h"
#include "CavityOnsiteOccupationsSquared.h"

std::vector<std::string> CavityLightIntensity::getHeader() const {
    return {"ada"};
}

std::vector<double> CavityLightIntensity::getValues() const {
    return {this->lightIntensity};
}

void CavityLightIntensity::calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObs) {
    const auto &onsiteOccupationsSquared
        = findObservable<CavityOnsiteOccupationsSquared>(primaryObs).getOccupationsSquared();
    this->lightIntensity = 0;
    for (std::size_t i{}; i < onsiteOccupationsSquared.size(); i++)
        for (std::size_t j{}; j < onsiteOccupationsSquared.size(); j++)
            this->lightIntensity += onsiteOccupationsSquared(i, j);
}