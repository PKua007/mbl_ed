//
// Created by Piotr Kubala on 29/01/2020.
//

#include <iterator>

#include "InverseParticipationRatio.h"

#include "utils/Quantity.h"
#include "simulation/RestorableHelper.h"


InverseParticipationRatio::InverseParticipationRatio(double relativeMiddleEnergy, double relativeMargin)
        : relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin}
{
    Expects(relativeMargin > 0);
    Expects(relativeMiddleEnergy - relativeMargin/2 > 0 && relativeMiddleEnergy + relativeMargin/2 < 1);
}

void InverseParticipationRatio::analyze(const Eigensystem &eigensystem, Logger &logger) {
    static_cast<void>(logger);

    Expects(eigensystem.hasEigenvectors());

    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    auto bandIndices = eigensystem.getIndicesOfNormalizedEnergiesInBand(this->relativeMiddleEnergy,
                                                                        this->relativeMargin);

    this->entries.reserve(bandIndices.size());
    for (auto i : bandIndices) {
        auto state = eigensystem.getEigenstate(i);
        double participationRatio = std::accumulate(state.begin(), state.end(), 0., [](double sum, double d) {
            return sum + d*d*d*d;
        });
        this->entries.emplace_back(normalizedEnergies[i], 1 / participationRatio);
    }
}

std::string InverseParticipationRatio::getName() const {
    return "ipr";
}

void InverseParticipationRatio::storeResult(std::ostream &out) const {
    std::copy(this->entries.begin(), this->entries.end(), std::ostream_iterator<Entry>(out, "\n"));
}

std::ostream &operator<<(std::ostream &out, const InverseParticipationRatio::Entry &entry) {
    out << entry.energy << " " << entry.ipr;
    return out;
}

void InverseParticipationRatio::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForVector(this->entries, binaryOut);
}

void InverseParticipationRatio::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForVector(this->entries, binaryIn);
}

void InverseParticipationRatio::clear() {
    this->entries.clear();
}
