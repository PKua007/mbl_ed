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
}

void EigenstateObservables::analyze(const Eigensystem &eigensystem, Logger &logger) {
    Expects(eigensystem.hasEigenvectors());

    std::size_t numBins = this->binEntries.size();
    for (std::size_t binIdx{}; binIdx < numBins; binIdx++) {
        double binBeg = static_cast<double>(binIdx) / numBins;
        double binEnd = static_cast<double>(binIdx + 1) / numBins;
        double binMid = (binBeg + binEnd) / 2;
        double binMargin = binEnd - binBeg;

        auto bandIndices = eigensystem.getIndicesOfNormalizedEnergiesInBand(binMid, binMargin);
        if (bandIndices.empty())
            continue;

        std::vector<double> observableValues(this->numValues);
        for (std::size_t bandIdx : bandIndices) {
            arma::vec doubleState = eigensystem.getEigenstate(bandIdx);
            arma::cx_vec state(arma::size(doubleState));
            std::copy(doubleState.begin(), doubleState.end(), state.begin());

            for (auto &primaryObservable : this->primaryObservables)
                primaryObservable->calculateForState(state);
            for (auto &secondaryObservable : this->secondaryObservables)
                secondaryObservable->calculateForObservables(this->primaryObservables);

            std::size_t offset{};
            for (const auto &storedObservable : this->storedObservables) {
                auto singleObservableValues = storedObservable->getValues();
                Assert(offset + singleObservableValues.size() <= this->numValues);
                std::transform(singleObservableValues.begin(), singleObservableValues.end(),
                               observableValues.begin() + offset, observableValues.begin() + offset,
                               std::plus<>{});
                offset += singleObservableValues.size();
            }
            Assert(offset == this->numValues);
        }

        for (std::size_t i{}; i < observableValues.size(); i++) {
            auto &binObservableValues = this->binEntries[binIdx].observableValues;
            Assert(binObservableValues.size() == observableValues.size());
            binObservableValues[i].push_back(observableValues[i]);
        }
    }
}

std::string EigenstateObservables::getName() const {
    return "obs";
}

void EigenstateObservables::storeResult(std::ostream &out) const {
    auto headerPlusError = [](const std::string &header) { return header + " d" + header; };
    out << "binStart ";
    std::transform(this->headers.begin(), this->headers.end(), std::ostream_iterator<std::string>(out, " "),
                   headerPlusError);
    out << std::endl;

    std::size_t numBins = this->binEntries.size();
    for (std::size_t binIdx{}; binIdx < numBins; binIdx++) {
        double binBeg = static_cast<double>(binIdx) / numBins;

        out << binBeg << " ";
        for (const auto &observableValues : this->binEntries[binIdx].observableValues) {
            Quantity binValue;
            binValue.calculateFromSamples(observableValues);
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

void EigenstateObservables::BinEntry::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForHistogram(this->observableValues, binaryOut);
}

void EigenstateObservables::BinEntry::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForHistogram(this->observableValues, binaryIn);
}

void EigenstateObservables::BinEntry::clear() {
    for (auto &observableValuesEntry : this->observableValues)
        observableValuesEntry.clear();
}
