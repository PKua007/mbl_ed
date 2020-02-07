//
// Created by pkua on 08.11.2019.
//

#include "MeanGapRatio.h"
#include "utils/Quantity.h"

MeanGapRatio::MeanGapRatio(double relativeMiddleEnergy, double relativeMargin)
        : relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin}
{
    Expects(relativeMargin > 0);
    Expects(relativeMiddleEnergy - relativeMargin/2 > 0 && relativeMiddleEnergy + relativeMargin/2 < 1);
}

void MeanGapRatio::analyze(const Eigensystem &eigensystem) {
    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    auto bandIndices = eigensystem.getIndicesOfNormalizedEnergiesInBand(this->relativeMiddleEnergy,
                                                                        this->relativeMargin);
    Assert(bandIndices.front() > 0);
    Assert(bandIndices.back() < normalizedEnergies.size() - 1);

    for (auto i : bandIndices) {
        double gap1 = normalizedEnergies[i] - normalizedEnergies[i - 1];
        double gap2 = normalizedEnergies[i + 1] - normalizedEnergies[i];
        this->gapRatios.push_back(gap1 < gap2 ? gap1/gap2 : gap2/gap1);
    }
}

Quantity MeanGapRatio::calculateMean() const {
    Quantity result;
    result.calculateFromSamples(this->gapRatios);
    return result;
}

std::string MeanGapRatio::getName() const {
    return "mgr";
}

std::vector<std::string> MeanGapRatio::getResultHeader() const {
    return {"meanGapRatio", "meanGapRatioError"};
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
