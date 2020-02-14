//
// Created by pkua on 01.11.2019.
//

#include "HamiltonianGenerator.h"
#include "utils/Assertions.h"

/**
 * @brief Returns vector after acting with b_{toSiteIndex}^\dagger b_{fromSiteIndex} on @a vector with a correct
 * constant
 */
std::optional<HopData>
HamiltonianGenerator::hoppingAction(const FockBase::Vector &fromVector, std::size_t fromSite,
                                    std::size_t toSite) const {
    Expects(toSite != fromSite);
    if (this->usePBC) {
        fromSite %= this->fockBase->getNumberOfSites();
        toSite %= this->fockBase->getNumberOfSites();
    } else if (fromSite >= this->fockBase->getNumberOfSites() || toSite >= this->fockBase->getNumberOfSites()) {
        return std::nullopt;
    }

    double constant = (fromVector[toSite] + 1) * fromVector[fromSite];
    if (constant == 0)
        return std::nullopt;

    HopData result;
    result.fromSite = fromSite;
    result.toSite = toSite;
    result.fromVector = fromVector;
    result.toVector = fromVector;
    result.toVector[fromSite]--;
    result.toVector[toSite]++;
    result.ladderConstant = std::sqrt(constant);

    return result;
}

arma::mat HamiltonianGenerator::generate() const {
    arma::mat result(this->fockBase->size(), this->fockBase->size(), arma::fill::zeros);

    for (std::size_t i = 0; i < this->fockBase->size(); i++) {
        for (auto &diagonalTerm : this->diagonalTerms)
            result(i, i) += diagonalTerm->calculate((*this->fockBase)[i], *this);

        for (std::size_t fromSite = 0; fromSite < this->fockBase->getNumberOfSites(); fromSite++) {
            auto hopData = this->hoppingAction((*this->fockBase)[i], fromSite, fromSite + 1);
            if (hopData == std::nullopt)
                continue;

            double hopConstant{};
            for (auto &hoppingTerm : this->hoppingTerms)
                hopConstant += hoppingTerm->calculate(*hopData, *this);

            hopConstant *= hopData->ladderConstant;
            std::size_t hoppedIndex = *(this->fockBase->findIndex(hopData->toVector));

            result(i, hoppedIndex) = hopConstant;
            result(hoppedIndex, i) = hopConstant;
        }

        for (std::size_t fromSite1 = 0; fromSite1 < this->fockBase->getNumberOfSites(); fromSite1++) {
            auto forward = this->hoppingAction((*this->fockBase)[i], fromSite1, fromSite1 + 1);
            auto backward = this->hoppingAction((*this->fockBase)[i], fromSite1 + 1, fromSite1);

            if (forward == std::nullopt && backward == std::nullopt)
                continue;

            for (std::size_t fromSite2 = 0; fromSite2 < this->fockBase->getNumberOfSites(); fromSite2++) {
                if (forward != std::nullopt) {
                    auto forwardForward = this->hoppingAction(forward->toVector, fromSite2, fromSite2 + 1);
                    if (forwardForward != std::nullopt) {
                        double hopConstant{};
                        for (auto &doubleHoppingTerm : this->doubleHoppingTerms)
                            hopConstant += doubleHoppingTerm->calculate(*forward, *forwardForward, *this);

                        hopConstant *= forward->ladderConstant * forwardForward->ladderConstant;
                        std::size_t hoppedIndex = *(this->fockBase->findIndex(forwardForward->toVector));

                        result(i, hoppedIndex) += hopConstant;
                        result(hoppedIndex, i) += hopConstant;
                    }

                    auto forwardBackward = this->hoppingAction(forward->toVector, fromSite2 + 1, fromSite2);
                    if (forwardBackward != std::nullopt) {
                        double hopConstant{};
                        for (auto &doubleHoppingTerm : this->doubleHoppingTerms)
                            hopConstant += doubleHoppingTerm->calculate(*forward, *forwardBackward, *this);

                        hopConstant *= forward->ladderConstant * forwardBackward->ladderConstant;
                        std::size_t hoppedIndex = *(this->fockBase->findIndex(forwardBackward->toVector));

                        result(hoppedIndex, i) += hopConstant;
                    }
                }

                if (backward != std::nullopt) {
                    auto backwardForward = this->hoppingAction(backward->toVector, fromSite2, fromSite2 + 1);
                    if (backwardForward != std::nullopt) {
                        double hopConstant{};
                        for (auto &doubleHoppingTerm : this->doubleHoppingTerms)
                            hopConstant += doubleHoppingTerm->calculate(*backward, *backwardForward, *this);

                        hopConstant *= backward->ladderConstant * backwardForward->ladderConstant;
                        std::size_t hoppedIndex = *(this->fockBase->findIndex(backwardForward->toVector));

                        result(hoppedIndex, i) += hopConstant;
                    }
                }
            }
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

std::vector<std::unique_ptr<DoubleHoppingTerm>> &HamiltonianGenerator::getDoubleHoppingTerms() {
    return this->doubleHoppingTerms;
}

void HamiltonianGenerator::addDiagonalTerm(std::unique_ptr<DiagonalTerm> term) {
    this->diagonalTerms.push_back(std::move(term));
}

void HamiltonianGenerator::addHoppingTerm(std::unique_ptr<HoppingTerm> term) {
    this->hoppingTerms.push_back(std::move(term));
}

void HamiltonianGenerator::addDoubleHoppingTerm(std::unique_ptr<DoubleHoppingTerm> term) {
    this->doubleHoppingTerms.push_back(std::move(term));
}

const FockBase &HamiltonianGenerator::getFockBase() const {
    return *(this->fockBase);
}
