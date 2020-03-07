//
// Created by Piotr Kubala on 19/02/2020.
//

#include "OccupationEvolution.h"

std::vector<OccupationEvolution::Occupations> OccupationEvolution::perform(double maxTime, std::size_t numSteps,
                                                                           std::size_t initialFockStateIdx,
                                                                           const Eigensystem &eigensystem,
                                                                           std::ostream &logger)
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
    logger << "[OccupationEvolution::perform] Calculating evolution operator... " << std::flush;
    arma::wall_clock timer;
    timer.tic();
    arma::cx_vec diagonalEvolution = arma::exp(-1i * dt * eigenenergies);
    arma::cx_mat fockBasisEvolution = eigenvectors * arma::diagmat(diagonalEvolution) * eigenvectors.t();
    logger << "done (" << timer.toc() << " s)." << std::endl;

    auto numOfParticlesObservables = prepareNumOfParticlesObservables(fockBase);
    auto numOfParticlesSquaresObservables = prepareNumOfParticlesSquaredObservables(fockBase);

    std::vector<Occupations> occupationEvolution = doPerformEvolution(numSteps, initialFockStateIdx, fockBase,
                                                                      fockBasisEvolution, numOfParticlesObservables,
                                                                      numOfParticlesSquaresObservables, logger);
    return occupationEvolution;
}

/**
 * @brief Based on the prepared evolution operator and observavles, do the actual evolution.
 */
std::vector<OccupationEvolution::Occupations>
OccupationEvolution::doPerformEvolution(std::size_t numSteps, size_t initialFockStateIdx, const FockBase &fockBase,
                                        const arma::cx_mat &fockBasisEvolution,
                                        const std::vector<arma::vec> &numOfParticlesObservables,
                                        const SymmetricMatrix<arma::vec> &numOfParticlesSquaredObservables,
                                        std::ostream &logger)
{
    arma::cx_vec evolvedState(fockBase.size(), arma::fill::zeros);
    evolvedState[initialFockStateIdx] = 1;

    arma::wall_clock timer;
    std::size_t numberOfSites = fockBase.getNumberOfSites();
    std::vector<Occupations> observablesEvolution = prepareOccupationVector(numSteps, numberOfSites);
    for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
        logger << "[OccupationEvolution::doPerformEvolution] Calculating expectation values for time step ";
        logger << timeIdx << "... " << std::flush;
        timer.tic();
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
        logger << "done (" << timer.toc() << " s). Evolving the state... " << std::flush;

        timer.tic();
        evolvedState = fockBasisEvolution * evolvedState;
        logger << "done (" << timer.toc() << " s)." << std::endl;
    }
    return observablesEvolution;
}

/**
 * @brief Prepare all pairs of n_i*n_j observables in the diagonal form and return SymmetricMatrix of them.
 */
SymmetricMatrix<arma::vec> OccupationEvolution::prepareNumOfParticlesSquaredObservables(const FockBase &fockBase) {
    std::size_t numberOfSites = fockBase.getNumberOfSites();
    SymmetricMatrix<arma::vec> result(numberOfSites);
    for (std::size_t site1Idx{}; site1Idx < numberOfSites; site1Idx++)
        for (std::size_t site2Idx = site1Idx; site2Idx < numberOfSites; site2Idx++)
            result(site1Idx, site2Idx) = numOfParticlesSquaredObservable(fockBase, site1Idx, site2Idx);
    return result;
}

/**
 * @brief Prepare all n_i observables in the diagonal form and return vector of them.
 */
std::vector<arma::vec> OccupationEvolution::prepareNumOfParticlesObservables(const FockBase &fockBase) {
    std::size_t numberOfSites = fockBase.getNumberOfSites();
    std::vector<arma::vec> result(numberOfSites);
    for (std::size_t siteIdx{}; siteIdx < numberOfSites; siteIdx++)
        result[siteIdx] = numOfParticlesObservable(fockBase, siteIdx);
    return result;
}

/**
 * @brief Initialize a vector of Occupations for given number of time steps and number of sites.
 */
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

/**
 * @param Return the diagonal of n_siteIdx observable in the diagonal basis.
 */
arma::vec OccupationEvolution::numOfParticlesObservable(const FockBase &fockBase, std::size_t siteIdx) {
    arma::vec result(fockBase.size());
    for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++)
        result[fockIdx] = fockBase[fockIdx][siteIdx];
    return result;
}

/**
 * @param Return the diagonal of n_site1Idx * n_site2Idx observable in the diagonal basis.
 */
arma::vec OccupationEvolution::numOfParticlesSquaredObservable(const FockBase &fockBase, std::size_t site1Idx,
                                                               std::size_t site2Idx)
{
    arma::vec result(fockBase.size());
    for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++)
        result[fockIdx] = fockBase[fockIdx][site1Idx] * fockBase[fockIdx][site2Idx];
    return result;
}