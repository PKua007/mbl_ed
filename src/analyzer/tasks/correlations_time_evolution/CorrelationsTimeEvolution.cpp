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
        auto observablesEvolution = OccupationEvolution::perform(this->times, initialIdx, eigensystem);

        Assert(observablesEvolution.size() == evolution.timeEntries.size());
        for (std::size_t i{}; i < observablesEvolution.size(); i++) {
            auto &timeEntry = evolution.timeEntries[i];
            auto &observables = observablesEvolution[i];
            timeEntry.addObservables(observables);
        }
    }

    this->meanEntries++;
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
            out << evolution.timeEntries[i].getValue();
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

std::string CorrelationsTimeEvolution::Correlations::getValue(std::size_t meanEntries) const {
    std::ostringstream out;
    out << (this->G / meanEntries);
    return out.str();
}

std::string CorrelationsTimeEvolution::OnsiteFluctuations::getValue(std::size_t meanEntries) const {
    std::ostringstream out;
    out << (this->rho / meanEntries);
    return out.str();
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

std::string CorrelationsTimeEvolution::TimeEntry::getValue() const {
    std::ostringstream out;
    out << this->t << " " << (this->x/this->meanEntries) << " ";

    auto valuePrinter = [this](const auto &entry) { return entry.getValue(meanEntries); };
    std::transform(this->correlations.begin(), this->correlations.end(), std::ostream_iterator<std::string>(out, " "),
                   valuePrinter);
    std::transform(this->borderlessCorrelations.begin(), this->borderlessCorrelations.end(),
                   std::ostream_iterator<std::string>(out, " "), valuePrinter);
    std::transform(this->onsiteFluctuations.begin(), this->onsiteFluctuations.end(),
                   std::ostream_iterator<std::string>(out, " "), valuePrinter);
    return out.str();
}

void CorrelationsTimeEvolution::TimeEntry::addObservables(const OccupationEvolution::Observables &observables) {
    this->meanEntries++;

    for (auto &onsiteFluctuation : this->onsiteFluctuations)
        onsiteFluctuation.addObservables(observables);
    for (auto &correlation : this->correlations)
        correlation.addObservables(observables, 0);
    for (auto &borderlessCorrelation : this->borderlessCorrelations)
        borderlessCorrelation.addObservables(observables, this->borderSize);

    this->x += 2*std::abs(std::accumulate(this->correlations.begin(), this->correlations.end(), 0.,
                                          [](double sum, const Correlations &corr) {
                                              return sum + corr.distance * corr.G;
                                          }));
}

CorrelationsTimeEvolution::TimeEntry::TimeEntry(double t, std::size_t borderSize, std::size_t numberOfSites)
        : t{t}, borderSize{borderSize}
{
    Expects(numberOfSites - 2*this->borderSize >= 2);

    this->onsiteFluctuations.resize(numberOfSites);
    std::size_t i{};
    for (auto &onsiteFluctuationsEntry : this->onsiteFluctuations)
        onsiteFluctuationsEntry.i = i++;

    this->correlations.resize(numberOfSites - 1);
    std::size_t d = 1;
    for (auto &correlationsEntry : this->correlations)
        correlationsEntry.distance = d++;

    this->borderlessCorrelations.resize(numberOfSites - 1 - 2 * this->borderSize);
    d = 1;
    for (auto &borderlessCorrelationEntry : this->borderlessCorrelations)
        borderlessCorrelationEntry.distance = d++;
}

std::size_t CorrelationsTimeEvolution::TimeEntry::getNumberOfSites() const {
    return 0;
}

std::string CorrelationsTimeEvolution::VectorEvolution::getHeader() const {
    Assert(!this->timeEntries.empty());
    std::ostringstream out;
    std::copy(this->initialVector.begin(), this->initialVector.end(), std::ostream_iterator<int>(out, "."));
    out.seekp(-1, std::ios_base::cur);
    out << "_" << this->timeEntries.front().getHeader();
    return out.str();
}