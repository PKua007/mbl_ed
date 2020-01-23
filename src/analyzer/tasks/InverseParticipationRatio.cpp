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

    double relativeFrom = this->relativeMiddleEnergy - this->relativeMargin/2;
    double relativeTo = this->relativeMiddleEnergy + this->relativeMargin/2;

    auto fromIt = std::lower_bound(normalizedEnergies.begin(), normalizedEnergies.end(), relativeFrom);
    auto toIt = std::lower_bound(normalizedEnergies.begin(), normalizedEnergies.end(), relativeTo);
    Assert(fromIt - normalizedEnergies.begin() >= 1);
    Assert(normalizedEnergies.end() - toIt >= 1);
    Assert(toIt - fromIt > 0);

    for (; fromIt < toIt; fromIt++) {
        auto vecIdx = fromIt - normalizedEnergies.begin();
        Assert(vecIdx >= 0 && static_cast<std::size_t>(vecIdx) < eigensystem.size());
        auto state = eigensystem.getEigenstate(vecIdx);
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
