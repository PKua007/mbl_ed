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
    if (this->getNumberOfSites() != 0)
        Expects(numberOfSites == this->getNumberOfSites());
    Expects(numberOfSites - 2*this->borderSize >= 2);

    if (this->getNumberOfSites() == 0) {
        Assert(!this->evolutions.empty());
        for (auto &evolution : this->evolutions) {
            Assert(!evolution.timeEntries.empty());
            for (auto &timeEntry : evolution.timeEntries) {
                timeEntry.onsiteFluctuations.resize(numberOfSites);
                std::size_t i{};
                for (auto &onsiteFluctuationsEntry : timeEntry.onsiteFluctuations)
                    onsiteFluctuationsEntry.i = i++;

                timeEntry.correlations.resize(numberOfSites - 1);
                std::size_t d = 1;
                for (auto &correlationsEntry : timeEntry.correlations)
                    correlationsEntry.distance = d++;

                timeEntry.borderlessCorrelations.resize(numberOfSites - 1 - 2 * this->borderSize);
                d = 1;
                for (auto &borderlessCorrelationEntry : timeEntry.borderlessCorrelations)
                    borderlessCorrelationEntry.distance = d++;
            }
        }
    }

    for (auto &evolution : this->evolutions) {
        std::size_t initialIdx = *fockBase.findIndex(evolution.initialVector);
        auto observablesEvolution = OccupationEvolution::perform(this->times, initialIdx, fockBase, eigensystem);

        Assert(observablesEvolution.size() == evolution.timeEntries.size());
        for (std::size_t i{}; i < observablesEvolution.size(); i++) {
            auto &timeEntry = evolution.timeEntries[i];
            auto &observables = observablesEvolution[i];

            for (auto &onsiteFluctuation : timeEntry.onsiteFluctuations)
                onsiteFluctuation.addObservables(observables);
            for (auto &correlations : timeEntry.correlations)
                correlations.addObservables(observables, 0);
            for (auto &borderlessCorrelations : timeEntry.borderlessCorrelations)
                borderlessCorrelations.addObservables(observables, this->borderSize);
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

    for (std::size_t i = 0; i < this->times.size(); i++) {
        for (const auto &evolution: this->evolutions)
            out << evolution.timeEntries[i];
        out << std::endl;
    }
}

CorrelationsTimeEvolution::CorrelationsTimeEvolution(double minTime, double maxTime, size_t numSteps,
                                                     CorrelationsTimeEvolution::TimeScaleType timeScaleType,
                                                     std::size_t borderSize,
                                                     const std::vector<FockBase::Vector> &vectorsToEvolve)
        : borderSize{borderSize}
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
    for (const auto &vectorToEvolve : vectorsToEvolve) {
        std::vector<TimeEntry> timeEntries;
        timeEntries.reserve(this->times.size());
        for (double time : this->times) {
            TimeEntry timeEntry;
            timeEntry.t = time;
            timeEntries.push_back(timeEntry);
        }
        this->evolutions.push_back({vectorToEvolve, timeEntries});
    }
}

std::size_t CorrelationsTimeEvolution::getNumberOfSites() const {
    Assert(!this->evolutions.empty());
    Assert(!this->evolutions.front().timeEntries.empty());
    return this->evolutions.front().timeEntries.front().onsiteFluctuations.size();
}

void CorrelationsTimeEvolution::Correlations::addObservables(const OccupationEvolution::Observables &observables,
                                                             std::size_t borderSize)
{
    std::size_t numberOfSites = observables.ns.size();
    Expects(this->distance < numberOfSites - 2*borderSize);

    double Gsum{};
    for (std::size_t i = borderSize; i < numberOfSites - this->distance - borderSize; i++)
        Gsum += observables.nns(i, i + this->distance) - observables.ns[i]*observables.ns[i + this->distance];
    std::size_t meanEntries = numberOfSites - this->distance - 2*borderSize;
    this->G += Gsum / meanEntries;
}

std::string CorrelationsTimeEvolution::Correlations::getHeader() const {
    return "G_" + std::to_string(this->distance);
}

std::string CorrelationsTimeEvolution::OnsiteFluctuations::getHeader() const {
    return "rho_" + std::to_string(this->i);
}

std::ostream &operator<<(std::ostream &out, const CorrelationsTimeEvolution::Correlations &corelations) {
    return out << corelations.distance;
}

std::ostream &operator<<(std::ostream &out, const CorrelationsTimeEvolution::OnsiteFluctuations &onsiteFluctuations) {
    return out << onsiteFluctuations.rho;
}

void CorrelationsTimeEvolution::OnsiteFluctuations::addObservables(const OccupationEvolution::Observables &observables)
{
    Expects(this->i < observables.ns.size());

    this->rho += observables.nns(this->i, this->i) - std::pow(observables.ns[this->i], 2);
}

std::string CorrelationsTimeEvolution::TimeEntry::getHeader() const {
    std::ostringstream out;

    out << "t x ";
    auto headerPrinter = [](const auto &entry) { return entry.getHeader(); };
    std::transform(this->correlations.begin(), this->correlations.end(), std::ostream_iterator<std::string>(out, " "),
                   headerPrinter);
    std::transform(this->borderlessCorrelations.begin(), this->borderlessCorrelations.end(),
                   std::ostream_iterator<std::string>(out, " "), headerPrinter);
    std::transform(this->onsiteFluctuations.begin(), this->onsiteFluctuations.end(),
                   std::ostream_iterator<std::string>(out, " "), headerPrinter);

    return out.str();
}

std::ostream &operator<<(std::ostream &out, const CorrelationsTimeEvolution::TimeEntry &timeEntry) {
    out << timeEntry.t << " " << timeEntry.x << " ";
    std::copy(timeEntry.correlations.begin(), timeEntry.correlations.end(),
              std::ostream_iterator<CorrelationsTimeEvolution::Correlations>(out, " "));
    std::copy(timeEntry.borderlessCorrelations.begin(), timeEntry.borderlessCorrelations.end(),
              std::ostream_iterator<CorrelationsTimeEvolution::Correlations>(out, " "));
    std::copy(timeEntry.onsiteFluctuations.begin(), timeEntry.onsiteFluctuations.end(),
              std::ostream_iterator<CorrelationsTimeEvolution::OnsiteFluctuations>(out, " "));
    return out;
}

std::string CorrelationsTimeEvolution::VectorEvolution::getHeader() const {
    Assert(!this->timeEntries.empty());
    std::ostringstream out;
    std::copy(this->initialVector.begin(), this->initialVector.end(), std::ostream_iterator<int>(out, "."));
    out.seekp(-1, std::ios_base::cur);
    out << "_" << this->timeEntries.front().getHeader();
    return out.str();
}