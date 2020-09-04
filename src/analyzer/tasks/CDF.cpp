//
// Created by Piotr Kubala on 17/01/2020.
//

#include <ostream>
#include <numeric>
#include <functional>

#include "core/Eigensystem.h"
#include "utils/Assertions.h"
#include "CDF.h"

void CDF::analyze(const Eigensystem &eigensystem, std::ostream &logger) {
    static_cast<void>(logger);

    for (auto &cdfEntry : this->cdfTable)
        cdfEntry.push_back(0);

    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    std::size_t steps = this->cdfTable.size();
    for (auto energy : normalizedEnergies) {
        auto it = this->cdfTable.begin();
        while(static_cast<double>(it - this->cdfTable.begin()) / static_cast<double>(steps - 1) < energy)
            it++;
        std::for_each(it, this->cdfTable.end(), [](auto &entry) { entry.back()++; });
    }

    for (auto &cdfEntry : this->cdfTable)
        cdfEntry.back() /= eigensystem.size();
}

std::string CDF::getName() const {
    return "cdf";
}

void CDF::storeResult(std::ostream &out) const {
    std::size_t steps = this->cdfTable.size();
    for (std::size_t i = 0; i < steps; i++) {
        out << (static_cast<double>(i) / static_cast<double>(steps - 1));
        out << " ";
        auto &entry = this->cdfTable[i];
        out << std::accumulate(entry.begin(), entry.end(), 0.) / entry.size();
        out << std::endl;
    }
}

CDF::CDF(std::size_t bins) : cdfTable(bins) {
    Expects(bins >= 2);
}

void CDF::storeState(std::ostream &binaryOut) const {
    std::size_t numBins = this->cdfTable.size();
    binaryOut.write(reinterpret_cast<const char*>(&numBins), sizeof(numBins));
    Assert(binaryOut.good());
    for (const auto &binEntries : this->cdfTable) {
        std::size_t numEntries = binEntries.size();
        binaryOut.write(reinterpret_cast<const char*>(&numEntries), sizeof(numEntries));
        binaryOut.write(reinterpret_cast<const char*>(binEntries.data()), sizeof(binEntries[0]) * numEntries);
    }
    Assert(binaryOut.good());
}

void CDF::joinRestoredState(std::istream &binaryIn) {
    std::size_t numBinsRestored{};
    binaryIn.read(reinterpret_cast<char*>(&numBinsRestored), sizeof(numBinsRestored));
    Assert(binaryIn.good());
    Assert(numBinsRestored == this->cdfTable.size());
    for (std::size_t i{}; i < this->cdfTable.size(); i++) {
        std::size_t numEntriesRestored{};
        binaryIn.read(reinterpret_cast<char*>(&numEntriesRestored), sizeof(numEntriesRestored));
        std::vector<double> entriesRestored(numEntriesRestored);
        binaryIn.read(reinterpret_cast<char*>(entriesRestored.data()), sizeof(entriesRestored[0])*numEntriesRestored);
        Assert(binaryIn.good());
        std::vector<double> &binEntries = this->cdfTable[i];
        binEntries.reserve(binEntries.size() + numEntriesRestored);
        binEntries.insert(binEntries.end(), entriesRestored.begin(), entriesRestored.end());
    }
}

void CDF::clear() {
    for (auto &binEntries : this->cdfTable)
        binEntries.clear();
}
