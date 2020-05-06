//
// Created by pkua on 02.05.2020.
//

#include <complex>

#include "EDEvolver.h"
#include "utils/Assertions.h"

using namespace std::complex_literals;

void EDEvolver::prepareFor(const arma::sp_mat &hamiltonian, const arma::cx_vec &initialState, double tMax,
                           std::size_t steps) {
    Expects(tMax > 0);
    Expects(steps >= 2);

    this->t = 0;
    this->dt = tMax / static_cast<double>(steps - 1);
    this->currentState = initialState;

    arma::mat eigvec;
    arma::vec eigval;
    Assert(arma::eig_sym(eigval, eigvec, arma::mat(hamiltonian)));
    arma::cx_vec diagonalEvolution = arma::exp(-1i * tMax * eigval);
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
