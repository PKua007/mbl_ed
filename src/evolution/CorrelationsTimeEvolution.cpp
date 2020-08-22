//
// Created by pkua on 08.05.2020.
//

#include <iterator>
#include <utility>

#include "CorrelationsTimeEvolution.h"
#include "utils/Assertions.h"
#include "Evolver.h"

void CorrelationsTimeEvolution::addEvolution(Evolver &evolver, std::ostream &logger,
                                             const std::vector<arma::cx_vec> &externalVectors)
{
    logger << externalVectors.size() << std::endl;
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

        logger << "[CorrelationsTimeEvolution::addEvolution] Evolving vector " << evolution.getInitialVectorName();
        logger << std::endl;

        OccupationEvolution occupationEvolution(this->fockBase);
        auto observablesEvolution = occupationEvolution.perform(this->timeSegmentation, initialState, evolver, logger);

        Assert(observablesEvolution.size() == evolution.timeEntries.size());
        for (std::size_t i{}; i < observablesEvolution.size(); i++) {
            auto &timeEntry = evolution.timeEntries[i];
            auto &observables = observablesEvolution[i];
            timeEntry.addObservables(observables);
        }
    }

    Assert(externalVectorsCounter == externalVectors.size());
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

    std::ostringstream out;
    if (std::holds_alternative<FockBase::Vector>(initialVector))
        out << std::get<FockBase::Vector>(initialVector);
    else
        out << std::get<ExternalVector>(initialVector).name;

    out << "_" << this->timeEntries.front().getHeader();
    return out.str();
}

std::string CorrelationsTimeEvolution::VectorEvolution::getInitialVectorName() const {
    if (std::holds_alternative<FockBase::Vector>(this->initialVector)) {
        std::ostringstream nameStream;
        nameStream << std::get<FockBase::Vector>(this->initialVector);
        return nameStream.str();
    } else {    // holds alternative ExternalVector
        ExternalVector externalVector = std::get<ExternalVector>(this->initialVector);
        return externalVector.name;
    }
}
