//
// Created by pkua on 08.05.2020.
//

#include <iterator>
#include <utility>

#include "TimeEvolution.h"
#include "utils/Assertions.h"
#include "Evolver.h"
#include "simulation/RestorableHelper.h"

void TimeEvolution::addEvolution(Evolver &evolver, Logger &logger, const std::vector<arma::cx_vec> &externalVectors) {
    Expects(externalVectors.size() == this->countExternalVectors());

    std::size_t externalVectorsCounter{};
    for (auto &evolution : this->vectorEvolutions) {
        arma::cx_vec initialState;

        if (std::holds_alternative<FockBase::Vector>(evolution.initialVector)) {
            FockBase::Vector initialFockVector = std::get<FockBase::Vector>(evolution.initialVector);
            auto initialIdx = fockBase->findIndex(initialFockVector);
            Assert(initialIdx.has_value());
            initialState = arma::cx_vec(this->fockBase->size(), arma::fill::zeros);
            initialState[*initialIdx] = 1;
        } else {    // holds alternative ExternalVector
            Assert(externalVectorsCounter < externalVectors.size());
            initialState = externalVectors[externalVectorsCounter];
            externalVectorsCounter++;
        }

        logger.info() << "Evolving vector " << evolution.getInitialVectorName() << std::endl;
        auto observablesEvolution = this->occupationEvolution->perform(this->timeSegmentation, initialState, evolver,
                                                                       logger);
        Assert(observablesEvolution.size() == evolution.timeEntries.size());
        std::transform(evolution.timeEntries.begin(), evolution.timeEntries.end(), observablesEvolution.begin(),
                       evolution.timeEntries.begin(), std::plus{});
    }
}

void TimeEvolution::storeResult(std::ostream &out) const {
    Assert(!this->vectorEvolutions.empty());

    auto headerPrinter = [](const auto &entry) { return entry.getHeader(); };
    std::transform(this->vectorEvolutions.begin(), this->vectorEvolutions.end(),
                   std::ostream_iterator<std::string>(out, ""), headerPrinter);
    out << std::endl;

    std::size_t numSteps = this->vectorEvolutions.front().timeEntries.size();
    for (std::size_t timeIdx{}; timeIdx < numSteps; timeIdx++) {
        for (const auto &evolution: this->vectorEvolutions)
            out << evolution.timeEntries[timeIdx].toString();
        out << std::endl;
    }
}

TimeEvolution::TimeEvolution(const TimeEvolutionParameters &parameters,
                             std::unique_ptr<OservablesTimeEvolution> occupationEvolution)
        : fockBase{parameters.fockBase}, occupationEvolution{std::move(occupationEvolution)},
          timeSegmentation{parameters.timeSegmentation}
{
    Expects(!parameters.vectorsToEvolve.empty());

    std::size_t numOfStoredObservableValues = parameters.countStoredObservableValues();
    this->occupationEvolution->setPrimaryObservables(parameters.primaryObservables);
    this->occupationEvolution->setSecondaryObservables(parameters.secondaryObservables);
    this->occupationEvolution->setStoredObservables(parameters.storedObservables);

    this->vectorEvolutions.reserve(parameters.vectorsToEvolve.size());
    for (const auto &vectorToEvolve : parameters.vectorsToEvolve) {
        VectorEvolution evolution;
        evolution.initialVector = vectorToEvolve;
        evolution.observablesHeader = parameters.generateStoredObservablesHeader();

        double lastMaxTime{};
        // Prepare entries for each time segment
        for (const auto &timeSegment : this->timeSegmentation) {
            for (std::size_t timeIdx{}; timeIdx < timeSegment.numSteps; timeIdx++) {
                double time = lastMaxTime +
                              (timeSegment.maxTime - lastMaxTime) / static_cast<double>(timeSegment.numSteps) * timeIdx;
                evolution.timeEntries.emplace_back(time, numOfStoredObservableValues);
            }
            lastMaxTime = timeSegment.maxTime;
        }
        // The last step should be separately added
        evolution.timeEntries.emplace_back(lastMaxTime, numOfStoredObservableValues);

        this->vectorEvolutions.push_back(evolution);
    }
}

void TimeEvolution::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForStaticRestorableVector(this->vectorEvolutions, binaryOut);
}

void TimeEvolution::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForStaticRestorableVector(this->vectorEvolutions, binaryIn);
}

void TimeEvolution::clear() {
    for (auto &vectorEvolution : this->vectorEvolutions)
        vectorEvolution.clear();
}

std::size_t TimeEvolution::countExternalVectors() const {
    std::size_t count{};
    using ExternalVector = TimeEvolutionParameters::ExternalVector;
    for (const auto &vectorEvolution : this->vectorEvolutions)
        if (std::holds_alternative<ExternalVector>(vectorEvolution.initialVector))
            count++;
    return count;
}

std::string TimeEvolution::VectorEvolution::getHeader() const {
    Assert(!this->timeEntries.empty());

    std::ostringstream out;
    if (std::holds_alternative<FockBase::Vector>(initialVector))
        out << std::get<FockBase::Vector>(initialVector);
    else    // holds alternative ExternalVector
        out << std::get<ExternalVector>(initialVector).name;
    out << "_t " << this->observablesHeader;

    return out.str();
}

std::string TimeEvolution::VectorEvolution::getInitialVectorName() const {
    if (std::holds_alternative<FockBase::Vector>(this->initialVector)) {
        std::ostringstream nameStream;
        nameStream << std::get<FockBase::Vector>(this->initialVector);
        return nameStream.str();
    } else {    // holds alternative ExternalVector
        ExternalVector externalVector = std::get<ExternalVector>(this->initialVector);
        return externalVector.name;
    }
}

void TimeEvolution::VectorEvolution::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForStaticRestorableVector(this->timeEntries, binaryOut);
}

void TimeEvolution::VectorEvolution::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForStaticRestorableVector(this->timeEntries, binaryIn);
}

void TimeEvolution::VectorEvolution::clear() {
    for (auto &timeEntry : this->timeEntries)
        timeEntry.clear();
}
