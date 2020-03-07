//
// Created by Piotr Kubala on 17/02/2020.
//

#include <iterator>

#include "CorrelationsTimeEvolution.h"
#include "utils/Assertions.h"

void CorrelationsTimeEvolution::analyze(const Eigensystem &eigensystem, std::ostream &logger) {
    Expects(eigensystem.hasEigenvectors());
    Expects(eigensystem.hasFockBase());

    const auto &fockBase = eigensystem.getFockBase();
    for (const auto &evolution : this->vectorEvolutions) {
        bool vectorIndexFound = fockBase.findIndex(evolution.initialVector).has_value();
        Expects(vectorIndexFound);
    }

    std::size_t numberOfSites = fockBase.getNumberOfSites();
    if (this->numberOfSitesDetermined())
        Expects(numberOfSites == this->getNumberOfSites());
    Expects(numberOfSites - 2*this->marginSize >= 2);

    if (!this->numberOfSitesDetermined())
        this->prepareTimeEntriesForNumberOfSites(numberOfSites);

    for (auto &evolution : this->vectorEvolutions) {
        logger << "[CorrelationsTimeEvolution::analyze] Evolving vector ";
        logger << evolution.getInitialVectorSignature() << std::endl;
        std::size_t initialIdx = *fockBase.findIndex(evolution.initialVector);
        auto observablesEvolution = OccupationEvolution::perform(this->maxTime, this->numSteps, initialIdx,
                                                                 eigensystem, logger);

        Assert(observablesEvolution.size() == evolution.timeEntries.size());
        for (std::size_t i{}; i < observablesEvolution.size(); i++) {
            auto &timeEntry = evolution.timeEntries[i];
            auto &observables = observablesEvolution[i];
            timeEntry.addObservables(observables);
        }
    }
}

void CorrelationsTimeEvolution::prepareTimeEntriesForNumberOfSites(std::size_t numberOfSites) {
    Assert(!this->vectorEvolutions.empty());
    for (auto &evolution : vectorEvolutions) {
        for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
            double time = maxTime / (static_cast<double>(numSteps) - 1) * timeIdx;
            evolution.timeEntries.emplace_back(time, marginSize, numberOfSites);
        }
    }
}

std::string CorrelationsTimeEvolution::getName() const {
    return "evolution";
}

void CorrelationsTimeEvolution::storeResult(std::ostream &out) const {
    Assert(!this->vectorEvolutions.empty());

    auto headerPrinter = [](const auto &entry) { return entry.getHeader(); };
    std::transform(this->vectorEvolutions.begin(), this->vectorEvolutions.end(), std::ostream_iterator<std::string>(out, ""),
                   headerPrinter);
    out << std::endl;

    for (std::size_t timeIdx{}; timeIdx < this->numSteps; timeIdx++) {
        for (const auto &evolution: this->vectorEvolutions)
            out << evolution.timeEntries[timeIdx].toString();
        out << std::endl;
    }
}

CorrelationsTimeEvolution::CorrelationsTimeEvolution(double maxTime, std::size_t numSteps, std::size_t marginSize,
                                                     const std::vector<FockBase::Vector> &vectorsToEvolve)
        : marginSize{marginSize}, maxTime{maxTime}, numSteps{numSteps}
{
    Expects(maxTime > 0);
    Expects(numSteps >= 2);
    Expects(!vectorsToEvolve.empty());

    this->vectorEvolutions.reserve(vectorsToEvolve.size());
    // For now there are no time entries, because the number of sites will be known only after the first eigensystem
    // arrives
    for (const auto &vectorToEvolve : vectorsToEvolve)
        this->vectorEvolutions.push_back({vectorToEvolve, {}});
}

std::size_t CorrelationsTimeEvolution::getNumberOfSites() const {
    Assert(this->numberOfSitesDetermined());
    return this->vectorEvolutions.front().timeEntries.front().getNumberOfSites();
}

bool CorrelationsTimeEvolution::numberOfSitesDetermined() const {
    Assert(!this->vectorEvolutions.empty());
    return !this->vectorEvolutions.front().timeEntries.empty();
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
