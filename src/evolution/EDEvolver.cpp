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

    this->dt = tMax / static_cast<double>(steps - 1);
    this->currentState = initialState;
    this->step = 0;
    this->steps = steps;

    arma::mat eigvec = this->eigensystem.getEigenstates();
    arma::vec eigval = this->eigensystem.getEigenenergies();
    arma::cx_vec diagonalEvolution = arma::exp(-1i * this->dt * eigval);
    this->evolutionOperator = eigvec * arma::diagmat(diagonalEvolution) * eigvec.t();
}

void EDEvolver::evolve() {
    // Actually this->step == this->steps - 1 here will give 1 step too much, but do not throw for convenience of use
    Assert(this->step < this->steps);
    this->step++;
    this->currentState = this->evolutionOperator * this->currentState;
}

const arma::cx_vec &EDEvolver::getCurrentState() const {
    return this->currentState;
}

EDEvolver::EDEvolver(const Eigensystem &eigensystem) : eigensystem{eigensystem} {
    Expects(eigensystem.hasEigenvectors());
}
