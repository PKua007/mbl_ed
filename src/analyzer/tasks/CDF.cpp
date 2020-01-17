//
// Created by Piotr Kubala on 17/01/2020.
//

#include <algorithm>
#include <ostream>
#include <numeric>

#include "utils/Assertions.h"
#include "CDF.h"

std::vector<double> CDF::getNormalizedEigenenergies(const std::vector<double> &eigenenergies) const {
    double groundLevel = eigenenergies.front();
    double highestLevel = eigenenergies.back();
    Expects(groundLevel < highestLevel);

    std::vector<double> normalizedEnergies;
    normalizedEnergies.reserve(eigenenergies.size());
    auto normalizer = [groundLevel, highestLevel](auto energy) {
        return (energy - groundLevel) / (highestLevel - groundLevel);
    };
    std::transform(eigenenergies.begin(), eigenenergies.end(), std::back_inserter(normalizedEnergies), normalizer);
    return normalizedEnergies;
}

void CDF::analyze(const std::vector<double> &eigenenergies) {
    std::size_t steps = this->cdfTable.size();
    auto normalizedEnergies = this->getNormalizedEigenenergies(eigenenergies);

    using namespace std::placeholders;
    std::for_each(this->cdfTable.begin(), this->cdfTable.end(), [](auto &entry) { entry.push_back(0); });
    for (auto energy : normalizedEnergies) {
        auto it = this->cdfTable.begin();
        while(static_cast<double>(it - this->cdfTable.begin()) / static_cast<double>(steps - 1) < energy)
            it++;

        std::for_each(it, this->cdfTable.end(), [](auto &entry) { entry.back()++; });
    }
    std::for_each(this->cdfTable.begin(), this->cdfTable.end(), [&eigenenergies](auto &entry) {
        entry.back() /= eigenenergies.size();
    });
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
