//
// Created by Piotr Kubala on 31/07/2020.
//

#include "BulkMeanGapRatio.h"

#include "utils/Quantity.h"

void BulkMeanGapRatio::analyze(const Eigensystem &eigensystem, std::ostream &logger) {
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
    std::size_t numBins = this->gapRatios.size();
    Assert(numBins > 0);
    std::size_t numEntries = this->gapRatios.front().size();
    binaryOut.write(reinterpret_cast<const char*>(&numBins), sizeof(numBins));
    binaryOut.write(reinterpret_cast<const char*>(&numEntries), sizeof(numEntries));
    Assert(binaryOut.good());
    for (const auto &binValues : this->gapRatios)
        binaryOut.write(reinterpret_cast<const char*>(binValues.data()), sizeof(binValues[0]) * numEntries);
    Assert(binaryOut.good());
}

void BulkMeanGapRatio::joinRestoredState(std::istream &binaryIn) {
    std::size_t numBinsRestored{};
    std::size_t numEntriesRestored{};
    binaryIn.read(reinterpret_cast<char*>(&numBinsRestored), sizeof(numBinsRestored));
    binaryIn.read(reinterpret_cast<char*>(&numEntriesRestored), sizeof(numEntriesRestored));
    Assert(binaryIn.good());
    Assert(numBinsRestored == this->gapRatios.size());
    for (std::size_t i{}; i < this->gapRatios.size(); i++) {
        std::vector<double> entriesRestored(numEntriesRestored);
        binaryIn.read(reinterpret_cast<char*>(entriesRestored.data()), sizeof(entriesRestored[0])*numEntriesRestored);
        Assert(binaryIn.good());
        std::vector<double> &binEntries = this->gapRatios[i];
        binEntries.reserve(binEntries.size() + numEntriesRestored);
        binEntries.insert(binEntries.end(), entriesRestored.begin(), entriesRestored.end());
    }
}

void BulkMeanGapRatio::clear() {
    for (auto &binValues : this->gapRatios)
        binValues.clear();
}
