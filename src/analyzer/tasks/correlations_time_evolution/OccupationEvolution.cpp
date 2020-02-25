//
// Created by Piotr Kubala on 19/02/2020.
//

#include "OccupationEvolution.h"

std::vector<OccupationEvolution::Observables> OccupationEvolution::perform(double minTime, double maxTime,
                                                                        std::size_t numSteps,
                                                                           size_t initialIdx,
                                                                           const Eigensystem &eigensystem)
{
    Expects(eigensystem.hasFockBase());
    Expects(minTime < maxTime);
    Expects(numSteps >= 2);
#ifndef NDEBUG
    Expects(eigensystem.isOrthonormal());
#endif

    const auto &fockBase = eigensystem.getFockBase();
    const auto &eigenvectors = eigensystem.getEigenstates();
    const auto &eigenenergies = eigensystem.getEigenenergies();
    std::size_t numberOfSites = fockBase.getNumberOfSites();

    std::vector<Observables> observablesEvolution;
    observablesEvolution.resize(numSteps);
    for (auto &observables : observablesEvolution)
        observables = Observables(numberOfSites);

    double dt = (maxTime - minTime) / (numSteps - 1);
    using namespace std::complex_literals;

    arma::cx_vec diagonalEvolution = arma::exp(-1i * dt * eigenenergies);
    arma::cx_mat fockEvolution = eigenvectors * arma::diagmat(diagonalEvolution) * eigenvectors.t();

    std::vector<arma::vec> ns(numberOfSites);
    std::vector<std::vector<arma::vec>> nns(numberOfSites);
    for (auto &nnsRow : nns)
        nnsRow.resize(numberOfSites);

    for (std::size_t site{}; site < numberOfSites; site++)
        ns[site] = numberOfParticlesObservable(fockBase, eigenvectors, site);

    for (std::size_t site1{}; site1 < numberOfSites; site1++)
        for (std::size_t site2 = site1; site2 < numberOfSites; site2++)
            nns[site1][site2] = numberOfParticlesSquaredObservable(fockBase, eigenvectors, site1, site2);

    arma::cx_vec currentVector(fockBase.size(), arma::fill::zeros);
    currentVector[initialIdx] = 1;

    for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
        for (std::size_t site{}; site < numberOfSites; site++) {
            observablesEvolution[timeIdx].ns[site] = calculateObservableValue(ns[site], currentVector);
        }

        for (std::size_t site1{}; site1 < numberOfSites; site1++) {
            for (std::size_t site2 = site1; site2 < numberOfSites; site2++) {
                observablesEvolution[timeIdx].nns(site1, site2) = calculateObservableValue(nns[site1][site2],
                                                                                           currentVector);
            }
        }

        currentVector = fockEvolution * currentVector;
    }

    return observablesEvolution;
}

double OccupationEvolution::calculateObservableValue(const arma::vec &observable, const arma::cx_vec &state)
{
    arma::cx_double result = arma::as_scalar(state.t() * arma::diagmat(observable) * state);
    return result.real();
}

arma::vec OccupationEvolution::numberOfParticlesObservable(const FockBase &fockBase,
                                                                 const arma::mat &eigenvectors, std::size_t site)
{
    arma::vec vec(fockBase.size());
    for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++)
        vec[fockIdx] = fockBase[fockIdx][site];

    return vec;
}

arma::vec OccupationEvolution::numberOfParticlesSquaredObservable(const FockBase &fockBase,
                                                                        const arma::mat &eigenvectors,
                                                                        std::size_t site1, std::size_t site2)
{
    arma::vec vec(fockBase.size());
    for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++)
        vec[fockIdx] = fockBase[fockIdx][site1] * fockBase[fockIdx][site2];

    return vec;
}