//
// Created by pkua on 08.11.2019.
//

#include <algorithm>

#include "GapRatioCalculator.h"
#include "Quantity.h"

GapRatioCalculator::GapRatioCalculator(double relativeMiddleEnergy, double relativeMargin)
        : relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin}
{
    Expects(relativeMargin > 0);
    Expects(relativeMiddleEnergy - relativeMargin/2 > 0 && relativeMiddleEnergy + relativeMargin/2 < 1);
}

void GapRatioCalculator::addEigenenergies(const std::vector<double> &eigenenergies) {
    double groundLevel = eigenenergies.front();
    double highestLevel = eigenenergies.back();
    Expects(groundLevel < highestLevel);

    std::vector<double> normalizedEnergies;
    normalizedEnergies.reserve(eigenenergies.size());
    auto normalizer = [groundLevel, highestLevel](auto energy) {
        return (energy - groundLevel) / (highestLevel - groundLevel);
    };
    std::transform(eigenenergies.begin(), eigenenergies.end(), std::back_inserter(normalizedEnergies), normalizer);

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

Quantity GapRatioCalculator::calculateMean() {
    Quantity result;
    result.calculateFromSamples(this->gapRatios);
    this->gapRatios.clear();
    return result;
}
