//
// Created by pkua on 08.11.2019.
//

#include "MeanGapRatio.h"

#include "utils/Quantity.h"
#include "simulation/RestorableHelper.h"

void MeanGapRatio::analyze(const Eigensystem &eigensystem, Logger &logger) {
    static_cast<void>(logger);

    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    std::vector<std::size_t> bandIndices = this->extractor.getBandIndices(eigensystem, logger);
    if (!bandIndices.empty() && bandIndices.front() == 0)
        bandIndices.erase(bandIndices.begin());
    if (!bandIndices.empty() && bandIndices.back() == eigensystem.size() - 1)
        bandIndices.pop_back();
    Assert(!bandIndices.empty());

    double singleGapRatio{};
    std::size_t numEntries{};
    for (auto i : bandIndices) {
        double gap1 = normalizedEnergies[i] - normalizedEnergies[i - 1];
        double gap2 = normalizedEnergies[i + 1] - normalizedEnergies[i];
        singleGapRatio += (gap1 < gap2 ? gap1/gap2 : gap2/gap1);
        numEntries++;
    }
    if (numEntries > 0)
        this->gapRatios.push_back(singleGapRatio / numEntries);
}

Quantity MeanGapRatio::calculateMean() const {
    Quantity result;
    result.calculateFromSamples(this->gapRatios);
    return result;
}

std::string MeanGapRatio::getName() const {
    return "mgr";
}

std::vector<std::string> MeanGapRatio::getResultHeader() const {
    return {"meanGapRatio", "meanGapRatioError"};
}

std::vector<std::string> MeanGapRatio::getResultFields() const {
    Quantity result = this->calculateMean();
    result.separator = Quantity::Separator::SPACE;
    std::stringstream resultStream;
    resultStream << result;
    std::string quantityValue, quantityError;
    resultStream >> quantityValue >> quantityError;
    return {quantityValue, quantityError};
}

void MeanGapRatio::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForVector(this->gapRatios, binaryOut);
}

void MeanGapRatio::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForVector(this->gapRatios, binaryIn);
}

void MeanGapRatio::clear() {
    this->gapRatios.clear();
}