//
// Created by Piotr Kubala on 06/12/2020.
//

#include "CavityOnsiteOccupations.h"

void CavityOnsiteOccupations::calculateForState(const arma::cx_vec &state) {
    for (std::size_t siteIdx{}; siteIdx < this->numOfSites; siteIdx++) {
        this->n_i[siteIdx] = Observable::calculateExpectedValue(this->diagonalObservables[siteIdx], state)
                             * this->cavityLongInteractions->calculateCosineForSite(siteIdx);
    }
}

CavityOnsiteOccupations::CavityOnsiteOccupations(std::shared_ptr<FockBasis> fockBasis,
                                                 std::shared_ptr<CavityLongInteraction> cavityLongInteractions)
        : numOfSites{fockBasis->getNumberOfSites()}, diagonalObservables(numOfSites, arma::vec(fockBasis->size())),
          n_i(numOfSites), fockBasis{std::move(fockBasis)}, cavityLongInteractions{std::move(cavityLongInteractions)}
{
    for (std::size_t siteIdx{}; siteIdx < this->numOfSites; siteIdx++)
        for (std::size_t fockIdx{}; fockIdx < this->fockBasis->size(); fockIdx++)
            this->diagonalObservables[siteIdx][fockIdx] = (*this->fockBasis)[fockIdx][siteIdx];
}

std::vector<std::string> CavityOnsiteOccupations::getHeader() const {
    std::vector<std::string> headerStrings;
    for (std::size_t siteIdx = 0; siteIdx < numOfSites; siteIdx++)
        headerStrings.push_back("n_" + std::to_string(siteIdx + 1) + "_cos");
    return headerStrings;
}

std::vector<double> CavityOnsiteOccupations::getValues() const {
    return this->n_i;
}