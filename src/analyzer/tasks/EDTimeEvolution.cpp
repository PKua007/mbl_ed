//
// Created by Piotr Kubala on 17/02/2020.
//

#include "EDTimeEvolution.h"

#include "evolution/EDEvolver.h"
#include "core/Eigensystem.h"
#include "utils/Assertions.h"

void EDTimeEvolution::analyze(const Eigensystem &eigensystem, Logger &logger) {
    EDEvolver evolver(eigensystem);
    this->timeEvolution.addEvolution(evolver, logger);
}

std::string EDTimeEvolution::getName() const {
    return "evolution";
}

void EDTimeEvolution::storeResult(std::ostream &out) const {
    this->timeEvolution.storeResult(out);
}

EDTimeEvolution::EDTimeEvolution(const TimeEvolutionParameters &parameters,
                                 std::unique_ptr<OservablesTimeEvolution> observablesEvolution)
        : timeEvolution(parameters, std::move(observablesEvolution))
{ }

void EDTimeEvolution::storeState(std::ostream &binaryOut) const {
    this->timeEvolution.storeState(binaryOut);
}

void EDTimeEvolution::joinRestoredState(std::istream &binaryIn) {
    this->timeEvolution.joinRestoredState(binaryIn);
}

void EDTimeEvolution::clear() {
    this->timeEvolution.clear();
}
