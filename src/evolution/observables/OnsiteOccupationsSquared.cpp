//
// Created by Piotr Kubala on 21/09/2020.
//

#include "OnsiteOccupationsSquared.h"

OnsiteOccupationsSquared::OnsiteOccupationsSquared(std::size_t numOfSites, std::shared_ptr<FockBase> fockBase)
        : numOfSites{numOfSites}, diagonalObservables(numOfSites, arma::vec(fockBase->size())), n_iN_j(numOfSites),
          header(numOfSites), fockBase{std::move(fockBase)}
{
    Expects(numOfSites > 0);

    for (std::size_t i{}; i < this->numOfSites; i++) {
        for (std::size_t j = i; j < this->numOfSites; j++) {
            for (std::size_t fockIdx{}; fockIdx < this->fockBase->size(); fockIdx++) {
                this->diagonalObservables(i, j)[fockIdx] =
                    (*this->fockBase)[fockIdx][i] * (*this->fockBase)[fockIdx][j];
                this->header(i, j) = "n_" + std::to_string(i + 1) + "N_" + std::to_string(j + 1);
            }
        }
    }
}

std::vector<std::string> OnsiteOccupationsSquared::getHeader() const {
    return std::vector(this->header.begin(), this->header.end());
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
