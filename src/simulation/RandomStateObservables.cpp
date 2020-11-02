//
// Created by Piotr Kubala on 19/10/2020.
//


#include "RandomStateObservables.h"
#include "RestorableHelper.h"
#include "utils/Quantity.h"
#include "utils/Assertions.h"


void RandomStateObservables::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForHistogram(this->notNormalizedValues, binaryOut);
    RestorableHelper::storeStateForHistogram(this->normalizedValues, binaryOut);
}

void RandomStateObservables::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForHistogram(this->notNormalizedValues, binaryIn);
    RestorableHelper::joinRestoredStateForHistogram(this->normalizedValues, binaryIn);
}

void RandomStateObservables::clear() {
    for (auto &entry : this->notNormalizedValues)
        entry.clear();
    for (auto &entry : this->normalizedValues)
        entry.clear();
}

void RandomStateObservables::performSimulation([[maybe_unused]] std::size_t simulationIndex,
                                               [[maybe_unused]] std::size_t totalSimulations,
                                               Logger &logger)
{
    logger.verbose() << "Simulation " << simulationIndex << " started..." << std::endl;
    arma::wall_clock timer;
    timer.tic();

    std::size_t size = this->basis->size();
    arma::cx_vec state(size);
    std::generate(state.begin(), state.end(), [this, size]() { return this->nextGaussian(1./size); });

    this->addStateToObservables(state, this->notNormalizedValues);
    this->addStateToObservables(arma::normalise(state), this->normalizedValues);

    logger.info() << "Simulation " << simulationIndex << " done (" << timer.toc() << " s)." << std::endl;
}

void RandomStateObservables::addStateToObservables(const arma::cx_vec &state,
                                                   std::vector<std::vector<double>> &observables) const
{
    for (auto &primaryObservable : primaryObservables)
        primaryObservable->calculateForState(state);
    for (auto &secondaryObservable : secondaryObservables)
        secondaryObservable->calculateForObservables(primaryObservables);

    std::size_t offset{};
    for (const auto &storedObservable : storedObservables) {
        auto singleObservableValues = storedObservable->getValues();
        Assert(offset + singleObservableValues.size() <= numValues);
        for (std::size_t i{}; i < singleObservableValues.size(); i++)
            observables[i + offset].push_back(singleObservableValues[i]);
        offset += singleObservableValues.size();
    }
    Assert(offset == numValues);
}

std::vector<std::string> RandomStateObservables::getValues() const {
    std::vector<std::string> result;
    result.reserve(this->numValues * 4);
    for (const auto &observableValues : this->notNormalizedValues)
        this->appendObservable(result, observableValues);
    for (const auto &observableValues : this->normalizedValues)
        this->appendObservable(result, observableValues);
    return result;
}

void RandomStateObservables::appendObservable(std::vector<std::string> &result,
                                              const std::vector<double> &observableValues) const
{
    Quantity meanValue;
    meanValue.calculateFromSamples(observableValues);
    meanValue.separator = Quantity::SPACE;
    std::stringstream inOut;
    inOut << meanValue;
    std::string value, error;
    inOut >> value >> error;
    Assert(inOut);
    result.push_back(value);
    result.push_back(error);
}

double RandomStateObservables::nextGaussian(double var) {
    double u1, u2;
    do {
        u1 = this->rnd->getDouble();
        u2 = this->rnd->getDouble();
    } while (u1 <= std::numeric_limits<double>::epsilon());
    return std::sqrt(-2 * var * std::log(u1)) * cos(2*M_PI*u2);
}

std::size_t RandomStateObservables::countStoredObservableValues() const {
    std::size_t numOfValues{};
    for (const auto &storedObservable : this->storedObservables)
        numOfValues += storedObservable->getHeader().size();
    return numOfValues;
}

std::vector<std::string> RandomStateObservables::generateStoredObservablesHeader() const {
    std::vector<std::string> headerStrings;
    headerStrings.reserve(this->numValues * 4);
    for (const auto &storedObservable : this->storedObservables) {
        auto observableHeaderStrings = storedObservable->getHeader();
        for (const auto &singleObservableString : observableHeaderStrings) {
            headerStrings.push_back(singleObservableString);
            headerStrings.push_back("d" + singleObservableString);
        }
    }

    size_t size = headerStrings.size();
    for (size_t i{}; i < size; ++i)
        headerStrings.push_back(headerStrings[i] + "_norm");

    return headerStrings;
}