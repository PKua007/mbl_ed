//
// Created by pkua on 08.05.2020.
//

#include <iterator>
#include <utility>

#include "CorrelationsTimeEvolution.h"
#include "utils/Assertions.h"
#include "Evolver.h"

void CorrelationsTimeEvolution::addEvolution(Evolver &evolver, std::ostream &logger) {
    for (auto &evolution : this->vectorEvolutions) {
        logger << "[EDCorrelationsTimeEvolution::analyze] Evolving vector ";
        logger << evolution.getInitialVectorSignature() << std::endl;

        auto initialIdx = fockBase->findIndex(evolution.initialVector);
        Assert(initialIdx.has_value());
        OccupationEvolution occupationEvolution(this->fockBase);
        auto observablesEvolution = occupationEvolution.perform(this->timeSegmentation, *initialIdx, evolver, logger);

        Assert(observablesEvolution.size() == evolution.timeEntries.size());
        for (std::size_t i{}; i < observablesEvolution.size(); i++) {
            auto &timeEntry = evolution.timeEntries[i];
            auto &observables = observablesEvolution[i];
            timeEntry.addObservables(observables);
        }
    }
}

void CorrelationsTimeEvolution::storeResult(std::ostream &out) const {
    Assert(!this->vectorEvolutions.empty());

    auto headerPrinter = [](const auto &entry) { return entry.getHeader(); };
    std::transform(this->vectorEvolutions.begin(), this->vectorEvolutions.end(), std::ostream_iterator<std::string>(out, ""),
                   headerPrinter);
    out << std::endl;

    std::size_t numSteps = this->vectorEvolutions.front().timeEntries.size();
    for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
        for (const auto &evolution: this->vectorEvolutions)
            out << evolution.timeEntries[timeIdx].toString();
        out << std::endl;
    }
}

CorrelationsTimeEvolution::CorrelationsTimeEvolution(const CorrelationsTimeEvolutionParameters &parameters)
        : fockBase{parameters.fockBase}, marginSize{parameters.marginSize},
          numberOfSites{parameters.numberOfSites}, timeSegmentation{parameters.timeSegmentation}
{
    Expects(!parameters.vectorsToEvolve.empty());

    this->vectorEvolutions.reserve(parameters.vectorsToEvolve.size());
    for (const auto &vectorToEvolve : parameters.vectorsToEvolve) {
        VectorEvolution evolution;
        evolution.initialVector = vectorToEvolve;

        double lastMaxTime{};
        // Prepare entries for each time segment
        for (const auto &timeSegment : this->timeSegmentation) {
            for (std::size_t timeIdx{}; timeIdx < timeSegment.numSteps; timeIdx++) {
                double time = lastMaxTime + (timeSegment.maxTime - lastMaxTime)
                              / static_cast<double>(timeSegment.numSteps) * timeIdx;
                evolution.timeEntries.emplace_back(time, marginSize, numberOfSites);
            }
            lastMaxTime = timeSegment.maxTime;
        }
        // The last step should be separately added
        evolution.timeEntries.emplace_back(lastMaxTime, marginSize, numberOfSites);

        this->vectorEvolutions.push_back(evolution);
    }
}

std::size_t CorrelationsTimeEvolution::getNumberOfSites() const {
    return this->numberOfSites;
}

std::string CorrelationsTimeEvolution::VectorEvolution::getHeader() const {
    Assert(!this->timeEntries.empty());
    std::string initialVectorSignature = this->getInitialVectorSignature();

    std::ostringstream out;
    out << initialVectorSignature << "_" << this->timeEntries.front().getHeader();
    return out.str();
}

std::string CorrelationsTimeEvolution::VectorEvolution::getInitialVectorSignature() const {
    std::ostringstream out;
    std::copy(this->initialVector.begin(), this->initialVector.end(), std::ostream_iterator<int>(out, "."));
    std::string outStr = out.str();
    // Remove unwanted '.' at the end
    outStr.pop_back();
    return outStr;
}
