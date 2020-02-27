//
// Created by Piotr Kubala on 19/02/2020.
//

#include "OccupationEvolution.h"

std::vector<OccupationEvolution::Occupations> OccupationEvolution::perform(double maxTime, std::size_t numSteps,
                                                                           std::size_t initialFockStateIdx,
                                                                           const Eigensystem &eigensystem)
{
    Expects(eigensystem.hasFockBase());
    Expects(maxTime > 0);
    Expects(numSteps >= 2);
#ifndef NDEBUG
    Expects(eigensystem.isOrthonormal());
#endif

    const auto &fockBase = eigensystem.getFockBase();
    const auto &eigenvectors = eigensystem.getEigenstates();
    const auto &eigenenergies = eigensystem.getEigenenergies();

    double dt = maxTime / (static_cast<double>(numSteps) - 1);

    using namespace std::complex_literals;
    arma::cx_vec diagonalEvolution = arma::exp(-1i * dt * eigenenergies);
    arma::cx_mat fockBasisEvolution = eigenvectors * arma::diagmat(diagonalEvolution) * eigenvectors.t();

    auto numOfParticlesObservables = prepareNumOfParticlesObservables(fockBase);
    auto numOfParticlesSquaresObservables = prepareNumOfParticlesSquaredObservable(fockBase);

    std::vector<Occupations> occupationEvolution = doPerformEvolution(numSteps, initialFockStateIdx, fockBase,
                                                                      fockBasisEvolution, numOfParticlesObservables,
                                                                      numOfParticlesSquaresObservables);
    return occupationEvolution;
}

std::vector<OccupationEvolution::Occupations>
OccupationEvolution::doPerformEvolution(std::size_t numSteps, size_t initialFockStateIdx, const FockBase &fockBase,
                                        const arma::cx_mat &fockBasisEvolution,
                                        const std::vector<arma::vec> &numOfParticlesObservables,
                                        const SymmetricMatrix<arma::vec> &numOfParticlesSquaredObservables)
{
    arma::cx_vec evolvedState(fockBase.size(), arma::fill::zeros);
    evolvedState[initialFockStateIdx] = 1;

    std::size_t numberOfSites = fockBase.getNumberOfSites();
    std::vector<Occupations> observablesEvolution = prepareOccupationVector(numSteps, numberOfSites);
    for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
        for (std::size_t site{}; site < numberOfSites; site++) {
            observablesEvolution[timeIdx].numParticles[site]
                = calculateObservableExpectedValue(numOfParticlesObservables[site], evolvedState);
        }

        for (std::size_t site1{}; site1 < numberOfSites; site1++) {
            for (std::size_t site2 = site1; site2 < numberOfSites; site2++) {
                observablesEvolution[timeIdx].numParticlesSquared(site1, site2)
                    = calculateObservableExpectedValue(numOfParticlesSquaredObservables(site1, site2), evolvedState);
            }
        }

        evolvedState = fockBasisEvolution * evolvedState;
    }
    return observablesEvolution;
}

SymmetricMatrix<arma::vec> OccupationEvolution::prepareNumOfParticlesSquaredObservable(const FockBase &fockBase) {
    std::size_t numberOfSites = fockBase.getNumberOfSites();
    SymmetricMatrix<arma::vec> result(numberOfSites);
    for (std::size_t site1Idx{}; site1Idx < numberOfSites; site1Idx++)
        for (std::size_t site2Idx = site1Idx; site2Idx < numberOfSites; site2Idx++)
            result(site1Idx, site2Idx) = numOfParticlesSquaredObservable(fockBase, site1Idx, site2Idx);
    return result;
}

std::vector<arma::vec> OccupationEvolution::prepareNumOfParticlesObservables(const FockBase &fockBase) {
    std::size_t numberOfSites = fockBase.getNumberOfSites();
    std::vector<arma::vec> result(numberOfSites);
    for (std::size_t siteIdx{}; siteIdx < numberOfSites; siteIdx++)
        result[siteIdx] = numOfParticlesObservable(fockBase, siteIdx);
    return result;
}

std::vector<OccupationEvolution::Occupations> OccupationEvolution::prepareOccupationVector(std::size_t numSteps,
                                                                                           std::size_t numberOfSites)
{
    std::vector<Occupations> result;
    result.resize(numSteps);
    for (auto &specificTimeOccupations : result)
        specificTimeOccupations = Occupations(numberOfSites);
    return result;
}

double OccupationEvolution::calculateObservableExpectedValue(const arma::vec &observable, const arma::cx_vec &state) {
    arma::cx_double result = arma::as_scalar(state.t() * arma::diagmat(observable) * state);
    return result.real();
}

arma::vec OccupationEvolution::numOfParticlesObservable(const FockBase &fockBase, std::size_t siteIdx) {
    arma::vec result(fockBase.size());
    for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++)
        result[fockIdx] = fockBase[fockIdx][siteIdx];
    return result;
}

arma::vec OccupationEvolution::numOfParticlesSquaredObservable(const FockBase &fockBase, std::size_t site1Idx,
                                                               std::size_t site2Idx)
{
    arma::vec result(fockBase.size());
    for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++)
        result[fockIdx] = fockBase[fockIdx][site1Idx] * fockBase[fockIdx][site2Idx];
    return result;
}