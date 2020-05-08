//
// Created by pkua on 02.05.2020.
//

#include <complex>

#include "EDEvolver.h"
#include "utils/Assertions.h"

using namespace std::complex_literals;

void EDEvolver::prepareFor(const arma::cx_vec &initialState, double tMax,
                           std::size_t steps) {
    Expects(tMax > 0);
    Expects(steps >= 2);
    Expects(initialState.size() == this->eigensystem.size());

    this->t = 0;
    this->dt = tMax / static_cast<double>(steps - 1);
    this->currentState = initialState;

    arma::mat eigvec = this->eigensystem.getEigenstates();
    arma::vec eigval = this->eigensystem.getEigenenergies();
    arma::cx_vec diagonalEvolution = arma::exp(-1i * this->dt * eigval);
    this->evolutionOperator = eigvec * arma::diagmat(diagonalEvolution) * eigvec.t();
}

void EDEvolver::evolve() {
    this->t += dt;
    this->currentState = this->evolutionOperator * this->currentState;
}

const arma::cx_vec &EDEvolver::getCurrentState() const {
    return this->currentState;
}

double EDEvolver::getCurrentTime() const {
    return this->t;
}

EDEvolver::EDEvolver(Eigensystem eigensystem) : eigensystem{std::move(eigensystem)} {
    Expects(eigensystem.hasEigenvectors());
}
