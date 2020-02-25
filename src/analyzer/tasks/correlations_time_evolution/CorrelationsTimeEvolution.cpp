//
// Created by Piotr Kubala on 17/02/2020.
//

#include <iterator>

#include "CorrelationsTimeEvolution.h"
#include "utils/Assertions.h"

void CorrelationsTimeEvolution::analyze(const Eigensystem &eigensystem) {
    Expects(eigensystem.hasEigenvectors());
    Expects(eigensystem.hasFockBase());

    const auto &fockBase = eigensystem.getFockBase();
    std::size_t numberOfSites = fockBase.getNumberOfSites();
    if (this->hasTimeEntries())
        Expects(numberOfSites == this->getNumberOfSites());
    Expects(numberOfSites - 2*this->borderSize >= 2);

    if (!this->hasTimeEntries()) {
        Assert(!this->evolutions.empty());
        for (auto &evolution : this->evolutions)
            for (double time : this->times)
                evolution.timeEntries.emplace_back(time, this->borderSize, numberOfSites);
    }

    for (auto &evolution : this->evolutions) {
        std::size_t initialIdx = *fockBase.findIndex(evolution.initialVector);
        auto observablesEvolution = OccupationEvolution::perform(this->minTime, this->maxTime, this->numSteps, initialIdx, eigensystem);

        Assert(observablesEvolution.size() == evolution.timeEntries.size());
        for (std::size_t i{}; i < observablesEvolution.size(); i++) {
            auto &timeEntry = evolution.timeEntries[i];
            auto &observables = observablesEvolution[i];
            timeEntry.addObservables(observables);
        }
    }
}

std::string CorrelationsTimeEvolution::getName() const {
    return "corr";
}

void CorrelationsTimeEvolution::storeResult(std::ostream &out) const {
    Assert(!this->evolutions.empty());

    auto headerPrinter = [](const auto &entry) { return entry.getHeader(); };
    std::transform(this->evolutions.begin(), this->evolutions.end(), std::ostream_iterator<std::string>(out, ""),
                   headerPrinter);
    out << std::endl;

    for (std::size_t i = 0; i < this->numSteps; i++) {
        for (const auto &evolution: this->evolutions)
            out << evolution.timeEntries[i].getValue();
        out << std::endl;
    }
}

CorrelationsTimeEvolution::CorrelationsTimeEvolution(double minTime, double maxTime, size_t numSteps,
                                                     CorrelationsTimeEvolution::TimeScaleType timeScaleType,
                                                     std::size_t borderSize,
                                                     const std::vector<FockBase::Vector> &vectorsToEvolve)
        : borderSize{borderSize}, minTime{minTime}, maxTime{maxTime}, numSteps{numSteps}
{
    Expects(minTime < maxTime);
    Expects(numSteps >= 2);
    Expects(!vectorsToEvolve.empty());

    this->times.reserve(numSteps);
    if (timeScaleType == Linear) {
        double factor = (maxTime - minTime) / (static_cast<double>(numSteps) - 1);
        for (std::size_t step{}; step < numSteps; step++)
            this->times.push_back(minTime + static_cast<double>(step) * factor);
    } else if (timeScaleType == Logarithmic) {
        Expects(minTime > 0);
        double logMinTime = std::log(minTime);
        double logMaxTime = std::log(maxTime);
        double factor = (logMaxTime - logMinTime) / (static_cast<double>(numSteps) - 1);
        for (std::size_t step{}; step < numSteps; step++)
            this->times.push_back(std::exp(logMinTime + static_cast<double>(step) * factor));
    } else {
        throw std::runtime_error("Unknown time scale type");
    }

    this->evolutions.reserve(vectorsToEvolve.size());
    for (const auto &vectorToEvolve : vectorsToEvolve)
        this->evolutions.push_back({vectorToEvolve, {}});
}

std::size_t CorrelationsTimeEvolution::getNumberOfSites() const {
    Assert(this->hasTimeEntries());
    return this->evolutions.front().timeEntries.front().getNumberOfSites();
}

bool CorrelationsTimeEvolution::hasTimeEntries() const {
    Assert(!this->evolutions.empty());
    return !this->evolutions.front().timeEntries.empty();
}

std::string CorrelationsTimeEvolution::VectorEvolution::getHeader() const {
    Assert(!this->timeEntries.empty());
    std::ostringstream out;
    std::copy(this->initialVector.begin(), this->initialVector.end(), std::ostream_iterator<int>(out, "."));
    out.seekp(-1, std::ios_base::cur);
    out << "_" << this->timeEntries.front().getHeader();
    return out.str();
}