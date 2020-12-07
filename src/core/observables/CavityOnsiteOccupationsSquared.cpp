//
// Created by Piotr Kubala on 06/12/2020.
//

#include "CavityOnsiteOccupationsSquared.h"

CavityOnsiteOccupationsSquared::
CavityOnsiteOccupationsSquared(std::shared_ptr<FockBasis> fockBasis,
                               std::shared_ptr<CavityLongInteraction> cavityLongInteractions)
        : numOfSites{fockBasis->getNumberOfSites()}, diagonalObservables(numOfSites, arma::vec(fockBasis->size())),
          n_iN_j(numOfSites), headerStrings(numOfSites), fockBasis{std::move(fockBasis)},
          cavityLongInteractions{std::move(cavityLongInteractions)}
{
    for (std::size_t i{}; i < this->numOfSites; i++) {
        for (std::size_t j = i; j < this->numOfSites; j++) {
            for (std::size_t fockIdx{}; fockIdx < this->fockBasis->size(); fockIdx++) {
                this->diagonalObservables(i, j)[fockIdx] =
                        (*this->fockBasis)[fockIdx][i] * (*this->fockBasis)[fockIdx][j];
                this->headerStrings(i, j) = "n_" + std::to_string(i + 1) + "N_" + std::to_string(j + 1) + "_cos";
            }
        }
    }
}

std::vector<std::string> CavityOnsiteOccupationsSquared::getHeader() const {
    return std::vector(this->headerStrings.begin(), this->headerStrings.end());
}

std::vector<double> CavityOnsiteOccupationsSquared::getValues() const {
    return std::vector(this->n_iN_j.begin(), this->n_iN_j.end());
}

void CavityOnsiteOccupationsSquared::calculateForState(const arma::cx_vec &state) {
    for (std::size_t i{}; i < this->numOfSites; i++)
        for (std::size_t j = i; j < this->numOfSites; j++)
            this->n_iN_j(i, j) = Observable::calculateExpectedValue(this->diagonalObservables(i, j), state)
                                 * this->cavityLongInteractions->calculateCosineForSite(i)
                                 * this->cavityLongInteractions->calculateCosineForSite(j);
}

const SymmetricMatrix<double> &CavityOnsiteOccupationsSquared::getOccupationsSquared() const {
    return this->n_iN_j;
}