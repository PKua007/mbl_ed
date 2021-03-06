//
// Created by Piotr Kubala on 17/01/2020.
//

#include <ostream>
#include <numeric>
#include <functional>

#include "core/Eigensystem.h"
#include "utils/Assertions.h"
#include "CDF.h"
#include "simulation/RestorableHelper.h"

void CDF::analyze(const Eigensystem &eigensystem, [[maybe_unused]] Logger &logger) {
    std::size_t numBins = this->cdfTable.size();
    std::vector<double> binsValue(numBins, 0);

    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    std::size_t steps = this->cdfTable.size();
    for (auto energy : normalizedEnergies) {
        auto it = binsValue.begin();
        while(static_cast<double>(it - binsValue.begin()) / static_cast<double>(steps - 1) < energy)
            it++;
        std::for_each(it, binsValue.end(), [](auto &binValue) { binValue++; });
    }

    for (std::size_t i{}; i < this->cdfTable.size(); i++)
        this->cdfTable[i].push_back(binsValue[i] / eigensystem.size());
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
    RestorableHelper::storeStateForHistogram(this->cdfTable, binaryOut);
}

void CDF::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForHistogram(this->cdfTable, binaryIn);
}

void CDF::clear() {
    for (auto &binEntries : this->cdfTable)
        binEntries.clear();
}
