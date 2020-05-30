//
// Created by Piotr Kubala on 19/02/2020.
//

#include "OccupationEvolution.h"

std::vector<OccupationEvolution::Occupations> OccupationEvolution::perform(const std::vector<EvolutionTimeSegment> &timeSegmentation,
                                                                           std::size_t initialFockStateIdx,
                                                                           const FockBase &fockBase, Evolver &evolver,
                                                                           std::ostream &logger)
{
    using namespace std::complex_literals;

    arma::cx_vec evolvedState(fockBase.size(), arma::fill::zeros);
    evolvedState[initialFockStateIdx] = 1;

    auto numOfParticlesObservables = prepareNumOfParticlesObservables(fockBase);
    auto numOfParticlesSquaresObservables = prepareNumOfParticlesSquaredObservables(fockBase);

    std::vector<Occupations> occupationEvolution;
    double lastMaxTime{};
    for (const auto &timeSegment : timeSegmentation) {
        logger << "[OccupationEvolution::perform] Calculating evolution operator... " << std::endl;
        arma::wall_clock timer;
        timer.tic();
        evolver.prepareFor(evolvedState, timeSegment.maxTime - lastMaxTime, timeSegment.numSteps + 1);
        logger << "[OccupationEvolution::perform] Calculating evolution operator done (" << timer.toc() << " s).";
        logger << std::endl;

        std::vector<Occupations> evolutionSegment = doPerformEvolution(timeSegment.numSteps, evolver, numOfParticlesObservables,
                                                                          numOfParticlesSquaresObservables, logger);
        occupationEvolution.insert(occupationEvolution.end(), evolutionSegment.begin(), evolutionSegment.end());
        evolvedState = evolver.getCurrentState();
        lastMaxTime = timeSegment.maxTime;
    }
    std::vector<Occupations> lastStep = doPerformEvolution(1, evolver, numOfParticlesObservables,
                                                                   numOfParticlesSquaresObservables, logger);
    occupationEvolution.push_back(lastStep.front());

    return occupationEvolution;
}

/**
 * @brief Based on the prepared evolution operator and observavles, do the actual evolution.
 */
std::vector<OccupationEvolution::Occupations>
OccupationEvolution::doPerformEvolution(std::size_t numSteps, Evolver &evolver,
                                        const std::vector<arma::vec> &numOfParticlesObservables,
                                        const SymmetricMatrix<arma::vec> &numOfParticlesSquaredObservables,
                                        std::ostream &logger)
{
    arma::wall_clock timer;
    std::size_t numberOfSites = numOfParticlesObservables.size();
    std::vector<Occupations> observablesEvolution = prepareOccupationVector(numSteps, numberOfSites);
    for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
        logger << "[OccupationEvolution::doPerformEvolution] Calculating expectation values for time step ";
        logger << timeIdx << "... " << std::flush;
        timer.tic();
        for (std::size_t site{}; site < numberOfSites; site++) {
            observablesEvolution[timeIdx].numParticles[site]
                = calculateObservableExpectedValue(numOfParticlesObservables[site], evolver.getCurrentState());
        }

        for (std::size_t site1{}; site1 < numberOfSites; site1++) {
            for (std::size_t site2 = site1; site2 < numberOfSites; site2++) {
                observablesEvolution[timeIdx].numParticlesSquared(site1, site2)
                    = calculateObservableExpectedValue(numOfParticlesSquaredObservables(site1, site2), evolver.getCurrentState());
            }
        }
        logger << "done (" << timer.toc() << " s). Evolving the state... " << std::flush;

        timer.tic();
        evolver.evolve();
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
    arma::vec modulusSquare(state.size());
    for (std::size_t i = 0; i < state.size(); i++)
        modulusSquare[i] = state[i].real() * state[i].real() + state[i].imag() * state[i].imag();

    return arma::as_scalar(modulusSquare.t() * observable);
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