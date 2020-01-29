//
// Created by Piotr Kubala on 22/01/2020.
//

#include "InverseParticipationRatio.h"
#include "utils/Quantity.h"

InverseParticipationRatio::InverseParticipationRatio(double relativeMiddleEnergy, double relativeMargin)
        : relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin}
{
    Expects(relativeMargin > 0);
    Expects(relativeMiddleEnergy - relativeMargin/2 > 0 && relativeMiddleEnergy + relativeMargin/2 < 1);
}

void InverseParticipationRatio::analyze(const Eigensystem &eigensystem) {
    Expects(eigensystem.hasEigenvectors());
    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    auto bandIndices = eigensystem.getIndicesOfNormalizedEnergiesInBand(this->relativeMiddleEnergy,
                                                                        this->relativeMargin);

    for (auto i : bandIndices) {
        auto state = eigensystem.getEigenstate(i);
        double participationRatio = std::accumulate(state.begin(), state.end(), 0., [](double sum, double d) {
            return sum + d*d*d*d;
        });
        this->ratios.push_back(1. / participationRatio);
    }
}

Quantity InverseParticipationRatio::calculateMean() const {
    Quantity result;
    result.calculateFromSamples(this->ratios);
    return result;
}

std::string InverseParticipationRatio::getName() const {
    return "ipr";
}

std::vector<std::string> InverseParticipationRatio::getResultHeader() const {
    return {"inverseParticipationRatio", "inverseParticipationRatioError"};
}

std::vector<std::string> InverseParticipationRatio::getResultFields() const {
    Quantity result = this->calculateMean();
    result.separator = Quantity::Separator::SPACE;
    std::stringstream resultStream;
    resultStream << result;
    std::string quantityValue, quantityError;
    resultStream >> quantityValue >> quantityError;
    return {quantityValue, quantityError};
}
