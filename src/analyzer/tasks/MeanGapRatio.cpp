//
// Created by pkua on 08.11.2019.
//

#include "MeanGapRatio.h"

#include "utils/Quantity.h"
#include "simulation/RestorableHelper.h"

MeanGapRatio::MeanGapRatio(double relativeMiddleEnergy, double relativeMargin)
        : relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin}
{
    Expects(relativeMargin > 0);
    Expects(relativeMiddleEnergy - relativeMargin/2 > 0 && relativeMiddleEnergy + relativeMargin/2 < 1);
}

MeanGapRatio::MeanGapRatio(const FockBasis::Vector &middleVector, double relativeMargin)
        : middleVector{middleVector}, relativeMargin{relativeMargin}
{
    Expects(relativeMargin > 0);
}

void MeanGapRatio::analyze(const Eigensystem &eigensystem, Logger &logger) {
    static_cast<void>(logger);

    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();

    double relativeMiddleEnergy_{};
    if (this->middleVector.has_value()) {
        relativeMiddleEnergy_ = calculateEnergyOfFockState(*this->middleVector, eigensystem);
        logger << "mgr calculated around: " << relativeMiddleEnergy_ << ". ";
    } else {
        relativeMiddleEnergy_ = this->relativeMiddleEnergy;
    }

    auto bandIndices = eigensystem.getIndicesOfNormalizedEnergiesInBand(relativeMiddleEnergy_, this->relativeMargin);
    if (!bandIndices.empty() && bandIndices.front() == 0)
        bandIndices.erase(bandIndices.begin());
    if (!bandIndices.empty() && bandIndices.back() == eigensystem.size() - 1)
        bandIndices.pop_back();
    Assert(!bandIndices.empty());

    double singleGapRatio{};
    std::size_t numEntries{};
    for (auto i : bandIndices) {
        double gap1 = normalizedEnergies[i] - normalizedEnergies[i - 1];
        double gap2 = normalizedEnergies[i + 1] - normalizedEnergies[i];
        singleGapRatio += (gap1 < gap2 ? gap1/gap2 : gap2/gap1);
        numEntries++;
    }
    if (numEntries > 0)
        this->gapRatios.push_back(singleGapRatio / numEntries);
}

double MeanGapRatio::calculateEnergyOfFockState(const FockBasis::Vector &state, const Eigensystem &eigensystem) const {
    Assert(eigensystem.hasFockBasis());
    Assert(eigensystem.hasEigenvectors());

    const auto &basis = eigensystem.getFockBasis();
    std::size_t basisIndex = *basis.findIndex(state);
    arma::vec fockVector(eigensystem.size(), arma::fill::zeros);
    fockVector[basisIndex] = 1;

    const arma::mat &eigvec = eigensystem.getEigenstates();
    const arma::vec &eigval = eigensystem.getEigenenergies();
    arma::vec diagonalVector = eigvec.t() * fockVector;

    double energy{};
    for (std::size_t i{}; i < diagonalVector.size(); i++)
        energy += std::pow(diagonalVector[i], 2) * eigval[i];

    double low = eigval.front();
    double high = eigval.back();
    double relativeEnergy = (energy - low) / (high - low);

    if (relativeEnergy - relativeMargin / 2 <= 0 || relativeEnergy + relativeMargin / 2 >= 1) {
        throw std::runtime_error("Margin " + std::to_string(this->relativeMargin) + " is to big around the given "
                                 "vector of the energy " + std::to_string(relativeEnergy));
    }
    return relativeEnergy;
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

void MeanGapRatio::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForVector(this->gapRatios, binaryOut);
}

void MeanGapRatio::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForVector(this->gapRatios, binaryIn);
}

void MeanGapRatio::clear() {
    this->gapRatios.clear();
}
