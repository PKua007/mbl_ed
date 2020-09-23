//
// Created by Piotr Kubala on 19/02/2020.
//

#include "OccupationEvolution.h"

#include <utility>

std::vector<CorrelationsTimeEntry>
OccupationEvolution::perform(const std::vector<EvolutionTimeSegment> &timeSegmentation,
                             const arma::cx_vec &initialState, Evolver &evolver, Logger &logger)
{
    this->timeStep = 0;
    this->time = 0;

    arma::cx_vec evolvedState = initialState;

    std::vector<CorrelationsTimeEntry> occupationEvolution;
    double lastMaxTime{};
    for (const auto &timeSegment : timeSegmentation) {
        logger.verbose() << "Calculating evolution operator... " << std::endl;
        arma::wall_clock timer;
        timer.tic();
        evolver.prepareFor(evolvedState, timeSegment.maxTime - lastMaxTime, timeSegment.numSteps);
        logger.info() << "Calculating evolution operator done (" << timer.toc() << " s).";
        logger << std::endl;

        std::vector<CorrelationsTimeEntry> evolutionSegment
            = this->performTimeSegmentEvolution(timeSegment.numSteps, evolver, logger);
        occupationEvolution.insert(occupationEvolution.end(), evolutionSegment.begin(), evolutionSegment.end());
        evolvedState = evolver.getCurrentState();
        lastMaxTime = timeSegment.maxTime;
    }
    std::vector<CorrelationsTimeEntry> lastEvolutionStep = performTimeSegmentEvolution(1, evolver, logger);
    occupationEvolution.push_back(lastEvolutionStep.front());

    return occupationEvolution;
}

/**
 * @brief Based on the prepared evolution operator and observavles, do the actual evolution of a single time segment
 * with a constant time step.
 */
std::vector<CorrelationsTimeEntry>
OccupationEvolution::performTimeSegmentEvolution(std::size_t numSteps, Evolver &evolver, Logger &logger) {
    arma::wall_clock timer;
    std::vector<CorrelationsTimeEntry> observablesEvolution;
    observablesEvolution.reserve(numSteps);
    for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
        logger.verbose() << "Calculating step " << this->timeStep << ", time " << this->time << " started...";
        logger << std::endl;

        timer.tic();

        for (auto &primaryObservable : this->primaryObservables)
            primaryObservable->calculateForState(evolver.getCurrentState());
        for (auto &secondaryObservable : this->secondaryObservables)
            secondaryObservable->calculateForObservables(this->primaryObservables);

        CorrelationsTimeEntry entry(this->time, this->numOfStoredValues);
        std::vector<double> values;
        values.reserve(this->numOfStoredValues);
        for (const auto &storedObservable : this->storedObservables) {
            auto observableValues = storedObservable->getValues();
            values.insert(values.end(), observableValues.begin(), observableValues.end());
        }
        entry.addValues(values);
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

OccupationEvolution::OccupationEvolution(std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables,
                                         std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables,
                                         std::vector<std::shared_ptr<Observable>> storedObservables)
        : primaryObservables{std::move(primaryObservables)}, secondaryObservables{std::move(secondaryObservables)},
          storedObservables{std::move(storedObservables)}
{
    for (const auto &storedObservable : this->storedObservables)
        this->numOfStoredValues += storedObservable->getHeader().size();
}