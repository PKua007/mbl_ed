//
// Created by pkua on 02.05.2020.
//

#include <complex>

#include "EDEvolver.h"
#include "utils/Assertions.h"

void EDEvolver::prepareFor(const arma::cx_vec &initialState, double maxTime, std::size_t numSteps_) {
    Expects(maxTime > 0);
    Expects(numSteps_ >= 2);
    Expects(initialState.size() == this->eigensystem.size());

    this->dt = maxTime / static_cast<double>(numSteps_ - 1);
    this->currentState = initialState;
    this->currentStep = 0;
    this->numSteps = numSteps_;

    using namespace std::complex_literals;
    arma::mat eigvec = this->eigensystem.getEigenstates();
    arma::vec eigval = this->eigensystem.getEigenenergies();
    arma::cx_vec diagonalEvolution = arma::exp(-1i * this->dt * eigval);
    this->evolutionOperator = eigvec * arma::diagmat(diagonalEvolution) * eigvec.t();
}

void EDEvolver::evolve() {
    // Actually this->currentStep == this->steps - 1 here will give 1 step too much, but do not throw for
    // convenience of use
    Assert(this->currentStep < this->numSteps);
    this->currentStep++;
    this->currentState = this->evolutionOperator * this->currentState;
}

const arma::cx_vec &EDEvolver::getCurrentState() const {
    return this->currentState;
}

EDEvolver::EDEvolver(const Eigensystem &eigensystem) : eigensystem{eigensystem} {
    Expects(eigensystem.hasEigenvectors());
}
