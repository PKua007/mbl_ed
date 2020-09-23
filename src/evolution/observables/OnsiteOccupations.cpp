//
// Created by Piotr Kubala on 21/09/2020.
//

#include "OnsiteOccupations.h"
#include "utils/Assertions.h"

void OnsiteOccupations::calculateForState(const arma::cx_vec &state) {
    for (std::size_t siteIdx{}; siteIdx < this->numOfSites; siteIdx++)
        this->n_i[siteIdx] = Observable::calculateExpectedValue(this->diagonalObservables[siteIdx], state);
}

OnsiteOccupations::OnsiteOccupations(std::shared_ptr<FockBase> fockBase)
        : numOfSites{fockBase->getNumberOfSites()}, diagonalObservables(numOfSites, arma::vec(fockBase->size())),
          n_i(numOfSites), fockBase{std::move(fockBase)}
{
    for (std::size_t siteIdx{}; siteIdx < this->numOfSites; siteIdx++)
        for (std::size_t fockIdx{}; fockIdx < this->fockBase->size(); fockIdx++)
            this->diagonalObservables[siteIdx][fockIdx] = (*this->fockBase)[fockIdx][siteIdx];
}

std::vector<std::string> OnsiteOccupations::getHeader() const {
    std::vector<std::string> headerStrings;
    for (std::size_t siteIdx = 0; siteIdx < numOfSites; siteIdx++)
        headerStrings.push_back("n_" + std::to_string(siteIdx + 1));
    return headerStrings;
}

std::vector<double> OnsiteOccupations::getValues() const {
    return this->n_i;
}
