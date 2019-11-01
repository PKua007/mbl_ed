//
// Created by pkua on 01.11.2019.
//

#include "HamiltonianGenerator.h"
#include "Assertions.h"

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
    arma::mat result(this->fockBase.size(), this->fockBase.size(), arma::fill::zeros);

    for (std::size_t i = 0; i < this->fockBase.size(); i++) {
        result(i, i) = this->getDiagonalElement(this->fockBase[i]);

        for (std::size_t j = 0; j < this->fockBase.getNumberOfSites() - 1; j++) {
            auto hoppedVector = this->hoppingAction(this->fockBase[i], j, j + 1);
            if (hoppedVector == std::nullopt)
                continue;

            auto [hoppedBase, hoppedConstant] = *hoppedVector;
            std::size_t hoppedIndex = *(this->fockBase.findIndex(hoppedBase));
            hoppedConstant *= this->getHoppingTerm(j, j + 1);

            result(i, hoppedIndex) = hoppedConstant;
            result(hoppedIndex, i) = hoppedConstant;
        }
    }

    return result;
}
