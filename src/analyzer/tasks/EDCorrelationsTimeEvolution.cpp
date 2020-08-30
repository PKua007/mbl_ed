//
// Created by Piotr Kubala on 17/02/2020.
//

#include "EDCorrelationsTimeEvolution.h"

#include "evolution/EDEvolver.h"
#include "core/Eigensystem.h"
#include "utils/Assertions.h"

void EDCorrelationsTimeEvolution::analyze(const Eigensystem &eigensystem, std::ostream &logger) {
    EDEvolver evolver(eigensystem);
    this->correlationsTimeEvolution.addEvolution(evolver, logger);
}

std::string EDCorrelationsTimeEvolution::getName() const {
    return "evolution";
}

void EDCorrelationsTimeEvolution::storeResult(std::ostream &out) const {
    this->correlationsTimeEvolution.storeResult(out);
}

EDCorrelationsTimeEvolution::EDCorrelationsTimeEvolution(const CorrelationsTimeEvolutionParameters &parameters)
        : correlationsTimeEvolution(parameters)
{ }