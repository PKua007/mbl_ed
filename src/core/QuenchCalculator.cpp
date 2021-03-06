//
// Created by Piotr Kubala on 18/08/2020.
//

#include "utils/Assertions.h"
#include "utils/Quantity.h"
#include "simulation/RestorableHelper.h"

#include "QuenchCalculator.h"

void QuenchCalculator::addQuench(const arma::sp_mat &initialHamiltonian, const arma::sp_mat &finalHamiltonian) {
    std::size_t numEigvals = std::min<std::size_t>(MIN_EIGVALS, initialHamiltonian.n_rows - 1);

    arma::vec initialMinEigvals;
    arma::mat initialMinEigvecs;
    Assert(arma::eigs_sym(initialMinEigvals, initialMinEigvecs, initialHamiltonian, numEigvals, "sa"));
    this->lastQuenchedState = initialMinEigvecs.col(0);

    arma::vec finalMinEigvals;
    arma::vec finalMaxEigvals;
    Assert(arma::eigs_sym(finalMinEigvals, finalHamiltonian, numEigvals, "sa"));
    Assert(arma::eigs_sym(finalMaxEigvals, finalHamiltonian, numEigvals, "la"));
    double Emin = finalMinEigvals.front();
    double Emax = finalMaxEigvals.back();

    arma::vec hamiltonianTimesQuenchedState = finalHamiltonian * this->lastQuenchedState;
    double quenchE = arma::as_scalar(this->lastQuenchedState.t() * hamiltonianTimesQuenchedState);
    double quenchE2 = arma::as_scalar(this->lastQuenchedState.t() * finalHamiltonian * hamiltonianTimesQuenchedState);
    double quenchEVariance = quenchE2 - quenchE * quenchE;
    // We just correct unmathematical negative values when it should be 0 originating from machine precision
    if (quenchEVariance < 0)
        quenchEVariance = 0;
    double quenchEpsilon = (quenchE - Emin) / (Emax - Emin);
    double quenchEpsilonVariance = quenchEVariance / std::pow((Emax - Emin), 2);

    this->quenchEpsilons.push_back(quenchEpsilon);
    this->quenchEpsilonVariances.push_back(quenchEpsilonVariance);
}

double QuenchCalculator::getMeanEpsilon() const {
    Quantity meanEpsilon;
    meanEpsilon.calculateFromSamples(this->quenchEpsilons);
    return meanEpsilon.value;
}

double QuenchCalculator::getMeanEpsilonQuantumUncertainty() const {
    Quantity meanVariance;
    meanVariance.calculateFromSamples(this->quenchEpsilonVariances);
    return std::sqrt(meanVariance.value);
}

double QuenchCalculator::getEpsilonAveragingSampleError() const {
    Quantity meanEpsilon;
    meanEpsilon.calculateFromSamples(this->quenchEpsilons);
    return meanEpsilon.error * std::sqrt(this->quenchEpsilons.size());
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

void QuenchCalculator::clear() {
    this->quenchEpsilons.clear();
    this->quenchEpsilonVariances.clear();
    this->lastQuenchedState.clear();
}

void QuenchCalculator::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForVector(this->quenchEpsilons, binaryOut);
    RestorableHelper::storeStateForVector(this->quenchEpsilonVariances, binaryOut);
    this->lastQuenchedState.save(binaryOut, arma::arma_binary);
    Assert(binaryOut);
}

void QuenchCalculator::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForVector(this->quenchEpsilons, binaryIn);
    RestorableHelper::joinRestoredStateForVector(this->quenchEpsilonVariances, binaryIn);
    this->lastQuenchedState.load(binaryIn, arma::arma_binary);
}
