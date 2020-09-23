//
// Created by Piotr Kubala on 17/02/2020.
//

#include "EDCorrelationsTimeEvolution.h"

#include "evolution/EDEvolver.h"
#include "core/Eigensystem.h"
#include "utils/Assertions.h"

void EDCorrelationsTimeEvolution::analyze(const Eigensystem &eigensystem, Logger &logger) {
    EDEvolver evolver(eigensystem);
    this->correlationsTimeEvolution.addEvolution(evolver, logger);
}

std::string EDCorrelationsTimeEvolution::getName() const {
    return "evolution";
}

void EDCorrelationsTimeEvolution::storeResult(std::ostream &out) const {
    this->correlationsTimeEvolution.storeResult(out);
}

EDCorrelationsTimeEvolution::EDCorrelationsTimeEvolution(const CorrelationsTimeEvolutionParameters &parameters,
                                                         std::unique_ptr<OccupationEvolution> occupationEvolution)
        : correlationsTimeEvolution(parameters, std::move(occupationEvolution))
{ }

void EDCorrelationsTimeEvolution::storeState(std::ostream &binaryOut) const {
    this->correlationsTimeEvolution.storeState(binaryOut);
}

void EDCorrelationsTimeEvolution::joinRestoredState(std::istream &binaryIn) {
    this->correlationsTimeEvolution.joinRestoredState(binaryIn);
}

void EDCorrelationsTimeEvolution::clear() {
    this->correlationsTimeEvolution.clear();
}
