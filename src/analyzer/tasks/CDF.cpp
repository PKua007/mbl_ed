//
// Created by Piotr Kubala on 17/01/2020.
//

#include <algorithm>
#include <ostream>
#include <numeric>
#include <functional>

#include "simulation/Eigensystem.h"
#include "utils/Assertions.h"
#include "CDF.h"

void CDF::analyze(const Eigensystem &eigensystem) {
    std::size_t steps = this->cdfTable.size();
    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();

    using namespace std::placeholders;
    std::for_each(this->cdfTable.begin(), this->cdfTable.end(), [](auto &entry) { entry.push_back(0); });
    for (auto energy : normalizedEnergies) {
        auto it = this->cdfTable.begin();
        while(static_cast<double>(it - this->cdfTable.begin()) / static_cast<double>(steps - 1) < energy)
            it++;

        std::for_each(it, this->cdfTable.end(), [](auto &entry) { entry.back()++; });
    }
    std::for_each(this->cdfTable.begin(), this->cdfTable.end(), [&eigensystem](auto &entry) {
        entry.back() /= eigensystem.size();
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
