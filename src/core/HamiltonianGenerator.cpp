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
HamiltonianGenerator::hoppingAction(const FockBasis::Vector &fromVector, std::size_t fromSite,
                                    std::size_t toSite) const {
    Expects(toSite != fromSite);
    if (this->usePBC) {
        fromSite %= this->fockBasis->getNumberOfSites();
        toSite %= this->fockBasis->getNumberOfSites();
    } else if (fromSite >= this->fockBasis->getNumberOfSites() || toSite >= this->fockBasis->getNumberOfSites()) {
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

arma::sp_mat HamiltonianGenerator::generate() const {
    arma::sp_mat result(this->fockBasis->size(), this->fockBasis->size());
    for (std::size_t vectorIdx = 0; vectorIdx < this->fockBasis->size(); vectorIdx++) {
        if (!this->diagonalTerms.empty())
            this->addDiagonalTerms(result, vectorIdx);
        for (const auto &hoppingTerm : this->hoppingTerms)
            this->addHoppingTerm(result, vectorIdx, *hoppingTerm);
        if (!this->doubleHoppingTerms.empty())
            this->addDoubleHoppingTerms(result, vectorIdx);
    }
    return result;
}

void HamiltonianGenerator::addDiagonalTerms(arma::sp_mat &result, std::size_t vectorIdx) const {
    for (auto &diagonalTerm : this->diagonalTerms)
        result(vectorIdx, vectorIdx) += diagonalTerm->calculate((*this->fockBasis)[vectorIdx], *this);
}

void HamiltonianGenerator::addHoppingTerm(arma::sp_mat &result, std::size_t fromIdx,
                                          const HoppingTerm &hoppingTerm) const
{
    for (std::size_t hoppingDistance : hoppingTerm.getHoppingDistances()) {
        for (std::size_t fromSite = 0; fromSite < this->fockBasis->getNumberOfSites(); fromSite++) {
            auto hopData = this->hoppingAction((*this->fockBasis)[fromIdx], fromSite, fromSite + hoppingDistance);
            if (hopData == std::nullopt)
                continue;

            double matrixElement{};
            matrixElement += hoppingTerm.calculate(*hopData, *this);

            matrixElement *= hopData->ladderConstant;
            std::size_t toIdx = *(this->fockBasis->findIndex(hopData->toVector));

            result(fromIdx, toIdx) += matrixElement;
            result(toIdx, fromIdx) += matrixElement;
        }
    }
}

auto HamiltonianGenerator::calculateDoubleHopMatrixElement(const HopData &firstHop, const HopData &secondHop) const {
    double matrixElement{};
    for (auto &doubleHoppingTerm : this->doubleHoppingTerms)
        matrixElement += doubleHoppingTerm->calculate(firstHop, secondHop, *this);

    matrixElement *= firstHop.ladderConstant * secondHop.ladderConstant;
    std::size_t toIdx = *(this->fockBasis->findIndex(secondHop.toVector));

    return std::make_tuple(matrixElement, toIdx);
}

void HamiltonianGenerator::performSecondHop(arma::sp_mat &result, std::size_t fromIdx, const HopData &firstHop) const {
    for (std::size_t fromSite2 = 0; fromSite2 < this->fockBasis->getNumberOfSites(); fromSite2++) {
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

void HamiltonianGenerator::addDoubleHoppingTerms(arma::sp_mat &result, std::size_t fromIdx) const {
    for (std::size_t fromSite1 = 0; fromSite1 < this->fockBasis->getNumberOfSites(); fromSite1++) {
        auto firstHopForward = this->hoppingAction((*this->fockBasis)[fromIdx], fromSite1, fromSite1 + 1);
        if (firstHopForward != std::nullopt)
            this->performSecondHop(result, fromIdx, *firstHopForward);

        auto firstHopBackward = this->hoppingAction((*this->fockBasis)[fromIdx], fromSite1 + 1, fromSite1);
        if (firstHopBackward != std::nullopt)
            this->performSecondHop(result, fromIdx, *firstHopBackward);
    }
}

size_t HamiltonianGenerator::getSiteDistance(std::size_t fromSite, std::size_t toIdx) const {
    Expects(fromSite < this->fockBasis->getNumberOfSites());
    Expects(toIdx < this->fockBasis->getNumberOfSites());
    std::size_t distance = abs(static_cast<int>(fromSite) - static_cast<int>(toIdx));
    if (this->usePBC && distance > this->fockBasis->getNumberOfSites() / 2)
        distance = this->fockBasis->getNumberOfSites() - distance;
    return distance;
}

const std::vector<std::shared_ptr<DiagonalTerm>> &HamiltonianGenerator::getDiagonalTerms() const {
    return this->diagonalTerms;
}

const std::vector<std::shared_ptr<HoppingTerm>> &HamiltonianGenerator::getHoppingTerms() const {
    return this->hoppingTerms;
}

const std::vector<std::shared_ptr<DoubleHoppingTerm>> &HamiltonianGenerator::getDoubleHoppingTerms() const {
    return this->doubleHoppingTerms;
}

void HamiltonianGenerator::addDiagonalTerm(std::shared_ptr<DiagonalTerm> term) {
    this->diagonalTerms.push_back(std::move(term));
}

void HamiltonianGenerator::addHoppingTerm(std::shared_ptr<HoppingTerm> term) {
    this->hoppingTerms.push_back(std::move(term));
}

void HamiltonianGenerator::addDoubleHoppingTerm(std::shared_ptr<DoubleHoppingTerm> term) {
    this->doubleHoppingTerms.push_back(std::move(term));
}

Eigensystem HamiltonianGenerator::calculateEigensystem(bool calculateEigenvectors) const {
    if (this->hoppingTerms.empty() && this->doubleHoppingTerms.empty()) {
        // For only diagonal terms there is no need to diagonalize
        arma::vec energies(this->fockBasis->size());
        for (auto &diagonalTerm : this->diagonalTerms)
            for (std::size_t i{}; i < this->fockBasis->size(); i++)
                energies[i] += diagonalTerm->calculate((*this->fockBasis)[i], *this);

        if (calculateEigenvectors)
            return Eigensystem(energies, arma::eye(this->fockBasis->size(), this->fockBasis->size()),
                               this->getFockBasis());
        else
            return Eigensystem(energies, this->getFockBasis());
    } else {
        // If off-diagonal elements are non-empty, diagonalization is needed
        arma::mat hamiltonian = arma::mat(this->generate());

        arma::vec armaEnergies;
        arma::mat armaEigvec;

        if (calculateEigenvectors) {
            Assert(arma::eig_sym(armaEnergies, armaEigvec, hamiltonian));
            return Eigensystem(armaEnergies, armaEigvec, this->getFockBasis());
        } else {
            Assert(arma::eig_sym(armaEnergies, hamiltonian));
            return Eigensystem(armaEnergies, this->getFockBasis());
        }
    }
}
