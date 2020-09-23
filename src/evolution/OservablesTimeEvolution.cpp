//
// Created by Piotr Kubala on 19/02/2020.
//

#include "OservablesTimeEvolution.h"

#include <utility>

std::vector<TimeEvolutionEntry>
OservablesTimeEvolution::perform(const std::vector<EvolutionTimeSegment> &timeSegmentation,
                                 const arma::cx_vec &initialState, Evolver &evolver, Logger &logger)
{
    this->timeStep = 0;
    this->time = 0;

    arma::cx_vec evolvedState = initialState;

    std::vector<TimeEvolutionEntry> observablesEvolution;
    double lastMaxTime{};
    for (const auto &timeSegment : timeSegmentation) {
        logger.verbose() << "Calculating evolution operator... " << std::endl;
        arma::wall_clock timer;
        timer.tic();
        evolver.prepareFor(evolvedState, timeSegment.maxTime - lastMaxTime, timeSegment.numSteps);
        logger.info() << "Calculating evolution operator done (" << timer.toc() << " s).";
        logger << std::endl;

        std::vector<TimeEvolutionEntry> evolutionSegment
            = this->performTimeSegmentEvolution(timeSegment.numSteps, evolver, logger);
        observablesEvolution.insert(observablesEvolution.end(), evolutionSegment.begin(), evolutionSegment.end());
        evolvedState = evolver.getCurrentState();
        lastMaxTime = timeSegment.maxTime;
    }
    std::vector<TimeEvolutionEntry> lastEvolutionStep = performTimeSegmentEvolution(1, evolver, logger);
    observablesEvolution.push_back(lastEvolutionStep.front());

    return observablesEvolution;
}

/**
 * @brief Based on the prepared evolution operator and observavles, do the actual evolution of a single time segment
 * with a constant time step.
 */
std::vector<TimeEvolutionEntry>
OservablesTimeEvolution::performTimeSegmentEvolution(std::size_t numSteps, Evolver &evolver, Logger &logger) {
    arma::wall_clock timer;
    std::vector<TimeEvolutionEntry> observablesEvolution;
    observablesEvolution.reserve(numSteps);
    for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
        logger.verbose() << "Calculating step " << this->timeStep << ", time " << this->time << " started...";
        logger << std::endl;

        timer.tic();

        for (auto &primaryObservable : this->primaryObservables)
            primaryObservable->calculateForState(evolver.getCurrentState());
        for (auto &secondaryObservable : this->secondaryObservables)
            secondaryObservable->calculateForObservables(this->primaryObservables);

        TimeEvolutionEntry entry(this->time, this->numOfObservableValues);
        std::vector<double> observableValues;
        observableValues.reserve(this->numOfObservableValues);
        for (const auto &storedObservable : this->storedObservables) {
            auto singleObservableValues = storedObservable->getValues();
            observableValues.insert(observableValues.end(), singleObservableValues.begin(),
                                    singleObservableValues.end());
        }
        entry.addValues(observableValues);
        observablesEvolution.push_back(entry);
        double observablesTime = timer.toc();

        timer.tic();
        evolver.evolve();
        double evolutionTime = timer.toc();

        logger.info() << "Calculating step " << this->timeStep << ", time " << this->time << " done (observables: ";
        logger << observablesTime << " s, state evolution: " << evolutionTime << " s)." << std::endl;

        this->timeStep++;
        this->time += evolver.getDt();
    }
    return observablesEvolution;
}

void OservablesTimeEvolution::setStoredObservables(const std::vector<std::shared_ptr<Observable>> &storedObservables_) {
    this->storedObservables = storedObservables_;
    this->numOfObservableValues = 0;
    for (const auto &storedObservable : this->storedObservables)
        this->numOfObservableValues += storedObservable->getHeader().size();
}
