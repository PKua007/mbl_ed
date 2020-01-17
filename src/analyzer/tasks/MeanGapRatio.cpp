//
// Created by pkua on 08.11.2019.
//

#include <algorithm>
#include <sstream>

#include "MeanGapRatio.h"
#include "utils/Quantity.h"

MeanGapRatio::MeanGapRatio(double relativeMiddleEnergy, double relativeMargin)
        : relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin}
{
    Expects(relativeMargin > 0);
    Expects(relativeMiddleEnergy - relativeMargin/2 > 0 && relativeMiddleEnergy + relativeMargin/2 < 1);
}

void MeanGapRatio::analyze(const std::vector<double> &eigenenergies) {
    std::vector<double> normalizedEnergies = this->getNormalizedEigenenergies(eigenenergies);

    double relativeFrom = this->relativeMiddleEnergy - this->relativeMargin/2;
    double relativeTo = this->relativeMiddleEnergy + this->relativeMargin/2;

    auto fromIt = std::lower_bound(normalizedEnergies.begin(), normalizedEnergies.end(), relativeFrom);
    auto toIt = std::lower_bound(normalizedEnergies.begin(), normalizedEnergies.end(), relativeTo);
    Assert(fromIt - normalizedEnergies.begin() >= 1);
    Assert(normalizedEnergies.end() - toIt >= 1);
    Assert(toIt - fromIt > 0);

    for (; fromIt < toIt; fromIt++) {
        double gap1 = *fromIt - *(fromIt - 1);
        double gap2 = *(fromIt + 1) - *fromIt;
        this->gapRatios.push_back(gap1 < gap2 ? gap1/gap2 : gap2/gap1);
    }
}

std::vector<double> MeanGapRatio::getNormalizedEigenenergies(const std::vector<double> &eigenenergies) const {
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

Quantity MeanGapRatio::calculateMean() const {
    Quantity result;
    result.calculateFromSamples(this->gapRatios);
    return result;
}

std::string MeanGapRatio::getName() const {
    return "mean_gap_ratio";
}

std::vector<std::string> MeanGapRatio::getResultHeader() const {
    return {"mean gap ratio", "mean gap ratio error"};
}

std::vector<std::string> MeanGapRatio::getResultFields() const {
    Quantity result = this->calculateMean();
    result.separator = Quantity::Separator::SPACE;
    std::stringstream resultStream;
    resultStream << result;
    std::string quantityValue, quantityError;
    resultStream >> quantityValue >> quantityError;
    return {quantityValue, quantityError};
}
