//
// Created by Piotr Kubala on 18/08/2020.
//

#include "utils/Assertions.h"
#include "utils/Quantity.h"

#include "QuenchCalculator.h"

void QuenchCalculator::addQuench(const arma::sp_mat &initialHamiltonian, const arma::sp_mat &finalHamiltonian) {
    std::size_t numEigvals = std::min<std::size_t>(MIN_EIGVALS, initialHamiltonian.n_rows - 1);

    arma::vec initialEigvals;
    arma::mat initialEigvecs;
    Assert(arma::eigs_sym(initialEigvals, initialEigvecs, initialHamiltonian, numEigvals, "sa"));
    this->lastQuenchedState = initialEigvecs.col(0);

    arma::vec finalMinimalEigvals;
    arma::vec finalMaximalEigvals;
    Assert(arma::eigs_sym(finalMinimalEigvals, finalHamiltonian, numEigvals, "sa"));
    Assert(arma::eigs_sym(finalMaximalEigvals, finalHamiltonian, numEigvals, "la"));
    double Emin = finalMinimalEigvals.front();
    double Emax = finalMaximalEigvals.back();

    arma::vec hamiltonianTimesQuenchedState = finalHamiltonian * this->lastQuenchedState;
    double Equench = arma::as_scalar(this->lastQuenchedState.t() * hamiltonianTimesQuenchedState);
    double E2quench = arma::as_scalar(this->lastQuenchedState.t() * finalHamiltonian * hamiltonianTimesQuenchedState);
    double varEquench = E2quench - Equench * Equench;
    double epsilonQuench = (Equench - Emin) / (Emax - Emin);
    double varEpsilonQuench = varEquench / std::pow((Emax - Emin), 2);

    this->quenchEpsilons.push_back(epsilonQuench);
    this->quenchEpsilonVariances.push_back(varEpsilonQuench);
}

double QuenchCalculator::getMeanEpsilon() const {
    Quantity quenchEpsilon;
    quenchEpsilon.calculateFromSamples(quenchEpsilons);
    return quenchEpsilon.value;
}

double QuenchCalculator::getMeanEpsilonQuantumUncertainty() const {
    Quantity quenchVariance;
    quenchVariance.calculateFromSamples(this->quenchEpsilonVariances);
    return std::sqrt(quenchVariance.value);
}

double QuenchCalculator::getEpsilonAveragingSampleError() const {
    Quantity quenchEpsilon;
    quenchEpsilon.calculateFromSamples(quenchEpsilons);
    return quenchEpsilon.error * std::sqrt(quenchEpsilons.size());
}

double QuenchCalculator::getLastQuenchEpsilon() const {
    if (this->quenchEpsilons.empty())
        return 0;
    else
        return this->quenchEpsilons.back();
}

double QuenchCalculator::getLastQuenchEpsilonQuantumUncertainty() const {
    if (this->quenchEpsilonVariances.empty())
        return 0;
    else
        return std::sqrt(this->quenchEpsilonVariances.back());
}
