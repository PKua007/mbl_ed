//
// Created by Piotr Kubala on 31/07/2020.
//

#include "BulkMeanGapRatio.h"

#include "simulation/RestorableHelper.h"
#include "utils/Quantity.h"

void BulkMeanGapRatio::analyze(const Eigensystem &eigensystem, Logger &logger) {
    static_cast<void>(logger);

    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();

    std::size_t numBins = this->gapRatios.size();
    for (std::size_t binIdx{}; binIdx < numBins; binIdx++) {
        double binBeg = static_cast<double>(binIdx) / numBins;
        double binEnd = static_cast<double>(binIdx + 1) / numBins;
        double binMid = (binBeg + binEnd) / 2;
        double binMargin = binEnd - binBeg;

        auto bandIndices = eigensystem.getIndicesOfNormalizedEnergiesInBand(binMid, binMargin);
        // Do not include ends of the spectrum - they do not have neighbours on both sides!
        if (!bandIndices.empty() && bandIndices.front() == 0)
            bandIndices.erase(bandIndices.begin());
        if (!bandIndices.empty() && bandIndices.back() == eigensystem.size() - 1)
            bandIndices.pop_back();

        if (bandIndices.empty())
            continue;

        for (auto i : bandIndices) {
            double gap1 = normalizedEnergies[i] - normalizedEnergies[i - 1];
            double gap2 = normalizedEnergies[i + 1] - normalizedEnergies[i];
            this->gapRatios[binIdx].push_back(gap1 < gap2 ? gap1 / gap2 : gap2 / gap1);
        }
    }
}

void BulkMeanGapRatio::storeResult(std::ostream &out) const {
    std::size_t numBins = this->gapRatios.size();
    for (std::size_t binIdx{}; binIdx < numBins; binIdx++) {
        double binBeg = static_cast<double>(binIdx) / numBins;
        Quantity binValue;
        binValue.calculateFromSamples(this->gapRatios[binIdx]);
        binValue.separator = Quantity::Separator::SPACE;
        out << binBeg << " " << binValue << std::endl;
    }
}

void BulkMeanGapRatio::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForHistogram(this->gapRatios, binaryOut);
}

void BulkMeanGapRatio::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForHistogram(this->gapRatios, binaryIn);
}

void BulkMeanGapRatio::clear() {
    for (auto &binValues : this->gapRatios)
        binValues.clear();
}
