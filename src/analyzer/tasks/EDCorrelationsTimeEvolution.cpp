//
// Created by Piotr Kubala on 17/02/2020.
//

#include "EDCorrelationsTimeEvolution.h"

#include "evolution/EDEvolver.h"
#include "simulation/Eigensystem.h"
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

EDCorrelationsTimeEvolution::EDCorrelationsTimeEvolution(double maxTime, std::size_t numSteps, std::size_t numberOfSites, std::size_t marginSize,
                                                         std::shared_ptr<FockBase> fockBase, const std::vector<FockBase::Vector> &vectorsToEvolve)
        : correlationsTimeEvolution(maxTime, numSteps, numberOfSites, marginSize, std::move(fockBase), vectorsToEvolve)
{ }