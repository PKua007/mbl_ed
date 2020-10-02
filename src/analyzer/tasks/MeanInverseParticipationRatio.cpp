//
// Created by Piotr Kubala on 22/01/2020.
//

#include "MeanInverseParticipationRatio.h"

#include "utils/Quantity.h"
#include "simulation/RestorableHelper.h"

MeanInverseParticipationRatio::MeanInverseParticipationRatio(double relativeMiddleEnergy, double relativeMargin)
        : relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin}
{
    Expects(relativeMargin > 0);
    Expects(relativeMiddleEnergy - relativeMargin/2 > 0 && relativeMiddleEnergy + relativeMargin/2 < 1);
}

void MeanInverseParticipationRatio::analyze(const Eigensystem &eigensystem, Logger &logger) {
    static_cast<void>(logger);

    Expects(eigensystem.hasEigenvectors());
    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    auto bandIndices = eigensystem.getIndicesOfNormalizedEnergiesInBand(this->relativeMiddleEnergy,
                                                                        this->relativeMargin);

    double singleRatio{};
    std::size_t numEntries{};
    for (auto i : bandIndices) {
        auto state = eigensystem.getEigenstate(i);
        double participationRatio = std::accumulate(state.begin(), state.end(), 0., [](double sum, double d) {
            return sum + d*d*d*d;
        });
        singleRatio += (1. / participationRatio);
        numEntries++;
    }
    if (numEntries > 0)
        this->ratios.push_back(singleRatio / numEntries);
}

Quantity MeanInverseParticipationRatio::calculateMean() const {
    Quantity result;
    result.calculateFromSamples(this->ratios);
    return result;
}

std::string MeanInverseParticipationRatio::getName() const {
    return "mipr";
}

std::vector<std::string> MeanInverseParticipationRatio::getResultHeader() const {
    return {"inverseParticipationRatio", "inverseParticipationRatioError"};
}

std::vector<std::string> MeanInverseParticipationRatio::getResultFields() const {
    Quantity result = this->calculateMean();
    result.separator = Quantity::Separator::SPACE;
    std::stringstream resultStream;
    resultStream << result;
    std::string quantityValue, quantityError;
    resultStream >> quantityValue >> quantityError;
    return {quantityValue, quantityError};
}

void MeanInverseParticipationRatio::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForVector(this->ratios, binaryOut);
}

void MeanInverseParticipationRatio::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForVector(this->ratios, binaryIn);
}

void MeanInverseParticipationRatio::clear() {
    this->ratios.clear();
}
