//
// Created by pkua on 01.11.2019.
//

#include "HamiltonianGenerator.h"
#include "utils/Assertions.h"

/**
 * @brief Returns vector after acting with b_{toSiteIndex}^\dagger b_{fromSiteIndex} on @a vector with a correct
 * constant
 */
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
        for (auto &diagonalTerm : this->diagonalTerms)
            result(i, i) += diagonalTerm->calculate((*this->fockBase)[i], *this);

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

            auto [hoppedBase, ladderConstant] = *hoppedVector;

            double hopConstant{};
            for (auto &hoppingTerm : this->hoppingTerms)
                hopConstant += hoppingTerm->calculate((*this->fockBase)[i], hoppedBase, from, to, *this);

            hopConstant *= ladderConstant;
            std::size_t hoppedIndex = *(this->fockBase->findIndex(hoppedBase));

            result(i, hoppedIndex) = hopConstant;
            result(hoppedIndex, i) = hopConstant;
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

std::vector<std::unique_ptr<DiagonalTerm>> &HamiltonianGenerator::getDiagonalTerms() {
    return this->diagonalTerms;
}

std::vector<std::unique_ptr<HoppingTerm>> &HamiltonianGenerator::getHoppingTerms() {
    return this->hoppingTerms;
}

void HamiltonianGenerator::addDiagonalTerm(std::unique_ptr<DiagonalTerm> term) {
    this->diagonalTerms.push_back(std::move(term));
}

void HamiltonianGenerator::addHoppingTerm(std::unique_ptr<HoppingTerm> term) {
    this->hoppingTerms.push_back(std::move(term));
}

const FockBase &HamiltonianGenerator::getFockBase() const {
    return *(this->fockBase);
}
