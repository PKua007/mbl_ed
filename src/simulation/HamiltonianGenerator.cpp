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

#include <iostream>
arma::mat HamiltonianGenerator::generate() const {
    arma::mat result(this->fockBase->size(), this->fockBase->size(), arma::fill::zeros);

    std::cout << "hamiltonian... " << std::flush;
    for (std::size_t vectorIdx = 0; vectorIdx < this->fockBase->size(); vectorIdx++) {
        this->addDiagonalTerms(result, vectorIdx);
        this->addHoppingTerms(result, vectorIdx);
        this->addDoubleHoppingTerms(result, vectorIdx);
    }
    std::cout << "done" << std::flush;

    return result;
}

void HamiltonianGenerator::addDiagonalTerms(arma::mat &result, std::size_t vectorIdx) const {
    for (auto &diagonalTerm : this->diagonalTerms)
        result(vectorIdx, vectorIdx) += diagonalTerm->calculate((*this->fockBase)[vectorIdx], *this);
}

void HamiltonianGenerator::addHoppingTerms(arma::mat &result, std::size_t fromIdx) const {
    for (std::size_t fromSite = 0; fromSite < this->fockBase->getNumberOfSites(); fromSite++) {
        auto hopData = this->hoppingAction((*this->fockBase)[fromIdx], fromSite, fromSite + 1);
        if (hopData == std::nullopt)
            continue;

        double matrixElement{};
        for (auto &hoppingTerm : this->hoppingTerms)
            matrixElement += hoppingTerm->calculate(*hopData, *this);

        matrixElement *= hopData->ladderConstant;
        std::size_t toIdx = *(this->fockBase->findIndex(hopData->toVector));

        result(fromIdx, toIdx) += matrixElement;
        result(toIdx, fromIdx) += matrixElement;
    }
}

auto HamiltonianGenerator::calculateDoubleHopMatrixElement(const HopData &firstHop, const HopData &secondHop) const {
    double matrixElement{};
    for (auto &doubleHoppingTerm : this->doubleHoppingTerms)
        matrixElement += doubleHoppingTerm->calculate(firstHop, secondHop, *this);

    matrixElement *= firstHop.ladderConstant * secondHop.ladderConstant;
    std::size_t toIdx = *(this->fockBase->findIndex(secondHop.toVector));

    return std::make_tuple(matrixElement, toIdx);
}

void HamiltonianGenerator::performSecondHop(arma::mat &result, std::size_t fromIdx, const HopData &firstHop) const {
    for (std::size_t fromSite2 = 0; fromSite2 < this->fockBase->getNumberOfSites(); fromSite2++) {
        auto secondHopForward = this->hoppingAction(firstHop.toVector, fromSite2, fromSite2 + 1);
        if (secondHopForward != std::nullopt) {
            auto[matrixElement, toIdx] = this->calculateDoubleHopMatrixElement(firstHop, *secondHopForward);
            result(toIdx, fromIdx) += matrixElement;
        }

        auto secondHopBackward = this->hoppingAction(firstHop.toVector, fromSite2 + 1, fromSite2);
        if (secondHopBackward != std::nullopt) {
            auto[matrixElement, toIdx] = this->calculateDoubleHopMatrixElement(firstHop, *secondHopBackward);
            result(toIdx, fromIdx) += matrixElement;
        }
    }
}

void HamiltonianGenerator::addDoubleHoppingTerms(arma::mat &result, std::size_t fromIdx) const {
    for (std::size_t fromSite1 = 0; fromSite1 < this->fockBase->getNumberOfSites(); fromSite1++) {
        auto firstHopForward = this->hoppingAction((*this->fockBase)[fromIdx], fromSite1, fromSite1 + 1);
        if (firstHopForward != std::nullopt)
            this->performSecondHop(result, fromIdx, *firstHopForward);

        auto firstHopBackward = this->hoppingAction((*this->fockBase)[fromIdx], fromSite1 + 1, fromSite1);
        if (firstHopBackward != std::nullopt)
            this->performSecondHop(result, fromIdx, *firstHopBackward);
    }
}

size_t HamiltonianGenerator::getSiteDistance(std::size_t fromSite, std::size_t toIdx) const {
    Expects(fromSite < this->fockBase->getNumberOfSites());
    Expects(toIdx < this->fockBase->getNumberOfSites());
    std::size_t distance = abs(static_cast<int>(fromSite) - static_cast<int>(toIdx));
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
