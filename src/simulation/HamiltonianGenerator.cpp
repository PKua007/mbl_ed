//
// Created by pkua on 01.11.2019.
//

#include "HamiltonianGenerator.h"
#include "utils/Assertions.h"

std::optional<std::pair<FockBase::Vector, double>>
HamiltonianGenerator::hoppingAction(const FockBase::Vector &vector, std::size_t fromSiteIndex,
                                    std::size_t toSiteIndex) const
{
    Expects(fromSiteIndex < vector.size());
    Expects(toSiteIndex < vector.size());

    double constant = (vector[toSiteIndex] + 1) * vector[fromSiteIndex];
    if (constant == 0)
        return std::nullopt;

    auto result = std::make_pair(vector, std::sqrt(constant));
    result.first[fromSiteIndex]--;
    result.first[toSiteIndex]++;
    return result;
}

arma::mat HamiltonianGenerator::generate() const {
    arma::mat result(this->fockBase->size(), this->fockBase->size(), arma::fill::zeros);

    for (std::size_t i = 0; i < this->fockBase->size(); i++) {
        result(i, i) = this->getDiagonalElement((*this->fockBase)[i]);

        for (std::size_t from = 0; from < this->fockBase->getNumberOfSites(); from++) {
            std::size_t to;
            if (this->usePBC)
                to = (from + 1) % this->fockBase->getNumberOfSites();
            else if (from == this->fockBase->getNumberOfSites() - 1)
                continue;
            else
                to = from + 1;

            auto hoppedVector = this->hoppingAction((*this->fockBase)[i], from, to);
            if (hoppedVector == std::nullopt)
                continue;

            auto [hoppedBase, hoppedConstant] = *hoppedVector;
            std::size_t hoppedIndex = *(this->fockBase->findIndex(hoppedBase));
            hoppedConstant *= this->getHoppingTerm(from, to);

            result(i, hoppedIndex) = hoppedConstant;
            result(hoppedIndex, i) = hoppedConstant;
        }
    }

    return result;
}

size_t HamiltonianGenerator::getSiteDistance(size_t fromSiteIndex, size_t toSiteIndex) const {
    Expects(fromSiteIndex < this->fockBase->getNumberOfSites());
    Expects(toSiteIndex < this->fockBase->getNumberOfSites());
    std::size_t distance = abs(static_cast<int>(fromSiteIndex) - static_cast<int>(toSiteIndex));
    if (this->usePBC && distance > this->fockBase->getNumberOfSites() / 2)
        distance = this->fockBase->getNumberOfSites() - distance;
    return distance;
}