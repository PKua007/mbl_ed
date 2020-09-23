//
// Created by Piotr Kubala on 21/09/2020.
//

#include "OnsiteOccupations.h"
#include "utils/Assertions.h"

void OnsiteOccupations::calculateForState(const arma::cx_vec &state) {
    for (std::size_t siteIdx{}; siteIdx < this->numOfSites; siteIdx++)
        this->n_i[siteIdx] = Observable::calculateExpectedValue(this->diagonalObservables[siteIdx], state);
}

OnsiteOccupations::OnsiteOccupations(std::size_t numOfSites, std::shared_ptr<FockBase> fockBase)
        : numOfSites{numOfSites}, diagonalObservables(numOfSites, arma::vec(fockBase->size())), n_i(numOfSites),
          fockBase{std::move(fockBase)}
{
    Expects(numOfSites > 0);

    for (std::size_t siteIdx{}; siteIdx < this->numOfSites; siteIdx++)
        for (std::size_t fockIdx{}; fockIdx < this->fockBase->size(); fockIdx++)
            this->diagonalObservables[siteIdx][fockIdx] = (*this->fockBase)[fockIdx][siteIdx];
}

std::vector<std::string> OnsiteOccupations::getHeader() const {
    std::vector<std::string> headerStrings;
    for (std::size_t i = 0; i < numOfSites; i++)
        headerStrings.push_back("n_" + std::to_string(i + 1));
    return headerStrings;
}

std::vector<double> OnsiteOccupations::getValues() const {
    return this->n_i;
}
