//
// Created by Piotr Kubala on 06/12/2020.
//

#include "CavityElectricField.h"
#include "CavityOnsiteOccupations.h"

std::vector<std::string> CavityElectricField::getHeader() const {
    return {"a"};
}

std::vector<double> CavityElectricField::getValues() const {
    return {this->electricField};
}

void CavityElectricField::calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObs) {
    auto onsiteOccupations = findObservable<CavityOnsiteOccupations>(primaryObs).getValues();
    this->electricField = std::accumulate(onsiteOccupations.begin(), onsiteOccupations.end(), 0.);
}
