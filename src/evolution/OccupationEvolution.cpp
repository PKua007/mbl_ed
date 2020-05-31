//
// Created by Piotr Kubala on 19/02/2020.
//

#include "OccupationEvolution.h"

std::vector<OccupationEvolution::Occupations>
OccupationEvolution::perform(const std::vector<EvolutionTimeSegment> &timeSegmentation, std::size_t initialFockStateIdx,
                             Evolver &evolver, std::ostream &logger)
{
    this->timeStep = 0;
    this->time = 0;

    arma::cx_vec evolvedState(this->fockBase->size(), arma::fill::zeros);
    evolvedState[initialFockStateIdx] = 1;

    std::vector<Occupations> occupationEvolution;
    double lastMaxTime{};
    for (const auto &timeSegment : timeSegmentation) {
        logger << "[OccupationEvolution::perform] Calculating evolution operator... " << std::endl;
        arma::wall_clock timer;
        timer.tic();
        evolver.prepareFor(evolvedState, timeSegment.maxTime - lastMaxTime, timeSegment.numSteps);
        logger << "[OccupationEvolution::perform] Calculating evolution operator done (" << timer.toc() << " s).";
        logger << std::endl;

        std::vector<Occupations> evolutionSegment = this->performTimeSegmentEvolution(timeSegment.numSteps, evolver,
                                                                                      logger);
        occupationEvolution.insert(occupationEvolution.end(), evolutionSegment.begin(), evolutionSegment.end());
        evolvedState = evolver.getCurrentState();
        lastMaxTime = timeSegment.maxTime;
    }
    std::vector<Occupations> lastEvolutionStep = performTimeSegmentEvolution(1, evolver, logger);
    occupationEvolution.push_back(lastEvolutionStep.front());

    return occupationEvolution;
}

/**
 * @brief Based on the prepared evolution operator and observavles, do the actual evolution of a single time segment
 * with a constant time step.
 */
std::vector<OccupationEvolution::Occupations>
OccupationEvolution::performTimeSegmentEvolution(std::size_t numSteps, Evolver &evolver, std::ostream &logger) {
    arma::wall_clock timer;
    std::size_t numberOfSites = this->fockBase->getNumberOfSites();
    std::vector<Occupations> observablesEvolution = this->prepareOccupationVector(numSteps, numberOfSites);
    for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
        logger << "[OccupationEvolution::performTimeSegmentEvolution] Calculating expectation values for step ";
        logger << this->timeStep << ", time " << this->time << "... " << std::flush;
        timer.tic();
        for (std::size_t site{}; site < numberOfSites; site++) {
            observablesEvolution[timeIdx].numParticles[site]
                = calculateObservableExpectedValue(this->numOfParticlesObservables[site], evolver.getCurrentState());
        }

        for (std::size_t site1{}; site1 < numberOfSites; site1++) {
            for (std::size_t site2 = site1; site2 < numberOfSites; site2++) {
                observablesEvolution[timeIdx].numParticlesSquared(site1, site2)
                    = calculateObservableExpectedValue(this->numOfParticlesSquaredObservables(site1, site2),
                                                       evolver.getCurrentState());
            }
        }
        logger << "done (" << timer.toc() << " s). Evolving the state... " << std::flush;

        timer.tic();
        evolver.evolve();
        logger << "done (" << timer.toc() << " s)." << std::endl;

        this->timeStep++;
        this->time += evolver.getDt();
    }
    return observablesEvolution;
}

/**
 * @brief Prepare all pairs of n_i*n_j observables in the diagonal form.
 */
void OccupationEvolution::prepareNumOfParticlesSquaredObservables() {
    std::size_t numberOfSites = this->fockBase->getNumberOfSites();
    this->numOfParticlesSquaredObservables = SymmetricMatrix<arma::vec>(numberOfSites);
    for (std::size_t site1Idx{}; site1Idx < numberOfSites; site1Idx++) {
        for (std::size_t site2Idx = site1Idx; site2Idx < numberOfSites; site2Idx++) {
            this->numOfParticlesSquaredObservables(site1Idx, site2Idx)
                = this->calculateNumOfParticlesSquaredObservable(site1Idx, site2Idx);
        }
    }
}

/**
 * @brief Prepare all n_i observables in the diagonal form.
 */
void OccupationEvolution::prepareNumOfParticlesObservables() {
    std::size_t numberOfSites = this->fockBase->getNumberOfSites();
    this->numOfParticlesObservables = std::vector<arma::vec>(numberOfSites);
    std::vector<arma::vec> result(numberOfSites);
    for (std::size_t siteIdx{}; siteIdx < numberOfSites; siteIdx++)
        this->numOfParticlesObservables[siteIdx] = this->calculateNumOfParticlesObservable(siteIdx);
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

double OccupationEvolution::calculateObservableExpectedValue(const arma::vec &observable,
                                                             const arma::cx_vec &state) const
{
    arma::vec modulusSquare(state.size());
    for (std::size_t i = 0; i < state.size(); i++)
        modulusSquare[i] = state[i].real() * state[i].real() + state[i].imag() * state[i].imag();

    return arma::as_scalar(modulusSquare.t() * observable);
}

/**
 * @param Return the diagonal of n_siteIdx observable in the diagonal basis.
 */
arma::vec OccupationEvolution::calculateNumOfParticlesObservable(std::size_t siteIdx) const {
    arma::vec result(this->fockBase->size());
    for (std::size_t fockIdx{}; fockIdx < this->fockBase->size(); fockIdx++)
        result[fockIdx] = (*this->fockBase)[fockIdx][siteIdx];
    return result;
}

/**
 * @param Return the diagonal of n_site1Idx * n_site2Idx observable in the diagonal basis.
 */
arma::vec OccupationEvolution::calculateNumOfParticlesSquaredObservable(std::size_t site1Idx,
                                                                        std::size_t site2Idx) const
{
    arma::vec result(this->fockBase->size());
    for (std::size_t fockIdx{}; fockIdx < this->fockBase->size(); fockIdx++)
        result[fockIdx] = (*this->fockBase)[fockIdx][site1Idx] * (*this->fockBase)[fockIdx][site2Idx];
    return result;
}

OccupationEvolution::OccupationEvolution(std::shared_ptr<FockBase> fockBase) : fockBase{std::move(fockBase)} {
    this->prepareNumOfParticlesObservables();
    this->prepareNumOfParticlesSquaredObservables();
}