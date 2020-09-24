//
// Created by Piotr Kubala on 21/09/2020.
//

#include "OnsiteOccupationsSquared.h"

OnsiteOccupationsSquared::OnsiteOccupationsSquared(std::shared_ptr<FockBase> fockBase)
        : numOfSites{fockBase->getNumberOfSites()}, diagonalObservables(numOfSites, arma::vec(fockBase->size())),
          n_iN_j(numOfSites), headerStrings(numOfSites), fockBase{std::move(fockBase)}
{
    for (std::size_t i{}; i < this->numOfSites; i++) {
        for (std::size_t j = i; j < this->numOfSites; j++) {
            for (std::size_t fockIdx{}; fockIdx < this->fockBase->size(); fockIdx++) {
                this->diagonalObservables(i, j)[fockIdx] =
                    (*this->fockBase)[fockIdx][i] * (*this->fockBase)[fockIdx][j];
                this->headerStrings(i, j) = "n_" + std::to_string(i + 1) + "N_" + std::to_string(j + 1);
            }
        }
    }
}

std::vector<std::string> OnsiteOccupationsSquared::getHeader() const {
    return std::vector(this->headerStrings.begin(), this->headerStrings.end());
}

std::vector<double> OnsiteOccupationsSquared::getValues() const {
    return std::vector(this->n_iN_j.begin(), this->n_iN_j.end());
}

void OnsiteOccupationsSquared::calculateForState(const arma::cx_vec &state) {
    for (std::size_t i{}; i < this->numOfSites; i++)
        for (std::size_t j = i; j < this->numOfSites; j++)
            this->n_iN_j(i, j) = Observable::calculateExpectedValue(this->diagonalObservables(i, j), state);
}

const SymmetricMatrix<double> &OnsiteOccupationsSquared::getOccupationsSquared() const {
    return this->n_iN_j;
}
