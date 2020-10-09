//
// Created by Piotr Kubala on 02/10/2020.
//

#include <algorithm>
#include <numeric>
#include <iterator>

#include "EigenstateObservables.h"

#include "simulation/RestorableHelper.h"
#include "utils/Assertions.h"
#include "utils/Quantity.h"

EigenstateObservables::EigenstateObservables(std::size_t numOfBins,
                                             std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables,
                                             std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables,
                                             std::vector<std::shared_ptr<Observable>> storedObservables)
        : binEntries(numOfBins), primaryObservables(std::move(primaryObservables)),
          secondaryObservables(std::move(secondaryObservables)), storedObservables(std::move(storedObservables))
{
    Expects(numOfBins > 0);
    this->numValues = this->countStoredObservableValues();
    for (auto &binEntry : this->binEntries)
        binEntry.observableValues.resize(this->numValues);
    this->header = "binStart " + this->generateStoredObservablesHeader();
}

auto EigenstateObservables::calculateBinRange(std::size_t binIdx, std::size_t numBins) {
    double binBeg = static_cast<double>(binIdx) / numBins;
    double binEnd = static_cast<double>(binIdx + 1) / numBins;

    // For first and last bin we want to go a little bit, respectively, lower that 0 and higher than 1 in order for
    // machine precision not to miss epsilon=0 and epsilon=1
    if (binIdx == 0)
        binBeg = -binEnd/2;
    if (binIdx == numBins - 1)
        binEnd = 1 + (1 - binBeg)/2;

    double binMid = (binBeg + binEnd) / 2;
    double binMargin = binEnd - binBeg;
    return std::make_pair(binMid, binMargin);
}

void EigenstateObservables::analyze(const Eigensystem &eigensystem, [[maybe_unused]] Logger &logger) {
    Expects(eigensystem.hasEigenvectors());

    arma::mat observables = this->calculateObservables(eigensystem);
    if (this->ostreamProvider != nullptr) {
        auto ostream = this->ostreamProvider->openOutputFile(
            this->fileSignature + "_" + std::to_string(this->eigensystemIdx) + "_obs.bin"
        );
        observables.save(*ostream, arma::arma_binary);
    }

    std::size_t numBins = this->binEntries.size();
    for (std::size_t binIdx{}; binIdx < numBins; binIdx++) {
        auto [binMid, binMargin] = this->calculateBinRange(binIdx, numBins);
        auto bandIndices = eigensystem.getIndicesOfNormalizedEnergiesInBand(binMid, binMargin);
        if (bandIndices.empty())
            continue;

        std::vector<double> observableValues = this->calculateMeanObservableValuesForBand(observables, bandIndices);
        for (std::size_t i{}; i < observableValues.size(); i++) {
            auto &binObservableValues = this->binEntries[binIdx].observableValues;
            Assert(binObservableValues.size() == observableValues.size());
            binObservableValues[i].push_back(observableValues[i]);
        }
    }

    this->eigensystemIdx++;
}

arma::mat EigenstateObservables::calculateObservables(const Eigensystem &eigensystem) const {
    arma::mat observables(eigensystem.size(), numValues);
    for (std::size_t i{}; i < eigensystem.size(); i++) {
        arma::vec doubleState = eigensystem.getEigenstate(i);
        arma::cx_vec state(arma::size(doubleState));
        std::copy(doubleState.begin(), doubleState.end(), state.begin());

        for (auto &primaryObservable : primaryObservables)
            primaryObservable->calculateForState(state);
        for (auto &secondaryObservable : secondaryObservables)
            secondaryObservable->calculateForObservables(primaryObservables);

        std::size_t offset{};
        for (const auto &storedObservable : storedObservables) {
            auto singleObservableValues = arma::rowvec(storedObservable->getValues());
            Assert(offset + singleObservableValues.size() <= numValues);
            observables(i, arma::span(offset, offset + singleObservableValues.size() - 1))
                = singleObservableValues;
            offset += singleObservableValues.size();
        }
        Assert(offset == numValues);
    }
    return observables;
}

std::vector<double>
EigenstateObservables::calculateMeanObservableValuesForBand(const arma::mat &observables,
                                                            const std::vector<std::size_t> &bandIndices) const
{
    std::vector<double> observableValues(numValues);
    for (std::size_t bandIdx : bandIndices) {
        std::transform(observables.begin_row(bandIdx), observables.end_row(bandIdx),
                       observableValues.begin(), observableValues.begin(),
                       std::plus<>{});
    }
    for (double &value : observableValues)
        value /= bandIndices.size();
    return observableValues;
}

std::string EigenstateObservables::getName() const {
    return "obs";
}

void EigenstateObservables::storeResult(std::ostream &out) const {
    out << this->header << std::endl;

    std::size_t numBins = this->binEntries.size();
    for (std::size_t binIdx{}; binIdx < numBins; binIdx++) {
        double binBeg = static_cast<double>(binIdx) / numBins;

        out << binBeg << " ";
        for (const auto &concreteObservableValues : this->binEntries[binIdx].observableValues) {
            Quantity binValue;
            binValue.calculateFromSamples(concreteObservableValues);
            binValue.separator = Quantity::Separator::SPACE;
            out << binValue << " ";
        }
        out << std::endl;
    }
}

void EigenstateObservables::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForStaticRestorableVector(this->binEntries, binaryOut);
}

void EigenstateObservables::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForStaticRestorableVector(this->binEntries, binaryIn);
}

void EigenstateObservables::clear() {
    for (auto &binEntry : this->binEntries)
        binEntry.clear();
}

std::size_t EigenstateObservables::countStoredObservableValues() const {
    std::size_t numOfValues{};
    for (const auto &storedObservable : this->storedObservables)
        numOfValues += storedObservable->getHeader().size();
    return numOfValues;
}

std::string EigenstateObservables::generateStoredObservablesHeader() const {
    std::vector<std::string> headerStrings;
    headerStrings.reserve(this->countStoredObservableValues());
    for (const auto &storedObservable : this->storedObservables) {
        auto observableHeaderStrings = storedObservable->getHeader();
        headerStrings.insert(headerStrings.end(), observableHeaderStrings.begin(), observableHeaderStrings.end());
    }

    std::ostringstream out;
    auto headerPlusError = [](const std::string &headerString) { return headerString + " d" + headerString; };
    std::transform(headerStrings.begin(), headerStrings.end(), std::ostream_iterator<std::string>(out, " "),
                   headerPlusError);
    return out.str();
}

void EigenstateObservables::BinEntry::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForHistogram(this->observableValues, binaryOut);
}

void EigenstateObservables::BinEntry::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForHistogram(this->observableValues, binaryIn);
}

void EigenstateObservables::BinEntry::clear() {
    for (auto &concreteObservableValues : this->observableValues)
        concreteObservableValues.clear();
}
