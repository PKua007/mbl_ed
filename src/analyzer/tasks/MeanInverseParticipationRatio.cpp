//
// Created by Piotr Kubala on 22/01/2020.
//

#include "MeanInverseParticipationRatio.h"
#include "utils/Quantity.h"

MeanInverseParticipationRatio::MeanInverseParticipationRatio(double relativeMiddleEnergy, double relativeMargin)
        : relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin}
{
    Expects(relativeMargin > 0);
    Expects(relativeMiddleEnergy - relativeMargin/2 > 0 && relativeMiddleEnergy + relativeMargin/2 < 1);
}

void MeanInverseParticipationRatio::analyze(const Eigensystem &eigensystem, std::ostream &logger) {
    static_cast<void>(logger);

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
    std::size_t ratiosSize = this->ratios.size();
    binaryOut.write(reinterpret_cast<const char*>(&ratiosSize), sizeof(ratiosSize));
    binaryOut.write(reinterpret_cast<const char*>(this->ratios.data()), sizeof(this->ratios[0]) * ratiosSize);
    Assert(binaryOut.good());
}

void MeanInverseParticipationRatio::joinRestoredState(std::istream &binaryIn) {
    std::size_t ratiosSizeRestored{};
    binaryIn.read(reinterpret_cast<char*>(&ratiosSizeRestored), sizeof(ratiosSizeRestored));
    Assert(binaryIn.good());
    std::vector<double> entriesRestored(ratiosSizeRestored);
    binaryIn.read(reinterpret_cast<char*>(entriesRestored.data()), sizeof(entriesRestored[0]) * ratiosSizeRestored);
    Assert(binaryIn.good());
    this->ratios.reserve(this->ratios.size() + ratiosSizeRestored);
    this->ratios.insert(this->ratios.end(), entriesRestored.begin(), entriesRestored.end());
}

void MeanInverseParticipationRatio::clear() {
    this->ratios.clear();
}
