//
// Created by Piotr Kubala on 01/11/2020.
//

#include <ostream>
#include <numeric>
#include <functional>

#include "core/Eigensystem.h"
#include "utils/Assertions.h"
#include "PDF.h"
#include "simulation/RestorableHelper.h"

void PDF::analyze(const Eigensystem &eigensystem, [[maybe_unused]] Logger &logger) {
    std::size_t numBins = this->pdfTable.size();
    std::vector<double> binsValue(numBins, 0);

    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    for (auto energy : normalizedEnergies) {
        auto binIdx = static_cast<std::size_t>(numBins * energy);
        Assert(binIdx <= numBins);
        // It is for energy = 1.0
        if (binIdx == numBins)
            binIdx = numBins - 1;
        binsValue[binIdx]++;
    }

    for (std::size_t i{}; i < this->pdfTable.size(); i++)
        this->pdfTable[i].push_back(binsValue[i] / eigensystem.size());
}

std::string PDF::getName() const {
    return "pdf";
}

void PDF::storeResult(std::ostream &out) const {
    std::size_t steps = this->pdfTable.size();
    for (std::size_t i = 0; i < steps; i++) {
        out << ((i + 0.5) / steps);
        out << " ";
        auto &entry = this->pdfTable[i];
        out << std::accumulate(entry.begin(), entry.end(), 0.) / entry.size();
        out << std::endl;
    }
}

PDF::PDF(std::size_t bins) : pdfTable(bins) {
    Expects(bins >= 2);
}

void PDF::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForHistogram(this->pdfTable, binaryOut);
}

void PDF::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForHistogram(this->pdfTable, binaryIn);
}

void PDF::clear() {
    for (auto &binEntries : this->pdfTable)
        binEntries.clear();
}

