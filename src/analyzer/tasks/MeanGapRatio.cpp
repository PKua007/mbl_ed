//
// Created by pkua on 08.11.2019.
//

#include "MeanGapRatio.h"

#include "utils/Quantity.h"
#include "simulation/RestorableHelper.h"

void MeanGapRatio::analyze(const Eigensystem &eigensystem, Logger &logger) {
    static_cast<void>(logger);

    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    std::vector<std::size_t> bandIndices = this->getBandIndices(eigensystem, normalizedEnergies, logger);
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

std::vector<std::size_t> MeanGapRatio::getBandIndices(const Eigensystem &eigensystem,
                                                      const arma::vec &normalizedEnergies, Logger &logger) const
{
    if (std::holds_alternative<VectorRange>(this->range)) {
        const auto &vectorRange = std::get<VectorRange>(this->range);
        double relativeMiddleEnergy = calculateEnergyOfFockState(vectorRange.middleVector, eigensystem);

        if (std::holds_alternative<WidthMargin>(vectorRange.epsilonMargin)) {
            double relativeMargin = std::get<WidthMargin>(vectorRange.epsilonMargin);
            if (relativeMiddleEnergy - relativeMargin / 2 <= 0 || relativeMiddleEnergy + relativeMargin / 2 >= 1) {
                throw std::runtime_error(
                    "Margin " + std::to_string(relativeMargin) + " is to big around the given vector of the energy "
                    + std::to_string(relativeMiddleEnergy)
                );
            }

            logger.info() << "Mgr calculated around: " << relativeMiddleEnergy << ". " << std::endl;
            return eigensystem.getIndicesOfNormalizedEnergiesInBand(relativeMiddleEnergy, relativeMargin);
        } else { // holds alternative NumberOfEnergiesMargin
            std::size_t numE = std::get<NumberOfEnergiesMargin>(vectorRange.epsilonMargin);
            auto bandIndices = eigensystem.getIndicesOfNumberOfNormalizedEnergies(relativeMiddleEnergy, numE);
            logger.info() << "Mgr calculated around: " << relativeMiddleEnergy << ", for: [";
            logger << normalizedEnergies[bandIndices.front()] << ", " << normalizedEnergies[bandIndices.back()] << "].";
            logger << std::endl;
            return bandIndices;
        }
    } else if (std::holds_alternative<EpsilonRange>(this->range)) {
        const auto &epsilonRange = std::get<EpsilonRange>(this->range);
        if (std::holds_alternative<WidthMargin>(epsilonRange.epsilonMargin)) {
            double margin = std::get<WidthMargin>(epsilonRange.epsilonMargin);
            return eigensystem.getIndicesOfNormalizedEnergiesInBand(epsilonRange.epsilonMiddle, margin);
        } else { // holds alternative NumberOfEnergiesMargin
            std::size_t numE = std::get<NumberOfEnergiesMargin>(epsilonRange.epsilonMargin);
            auto bandIndices = eigensystem.getIndicesOfNumberOfNormalizedEnergies(epsilonRange.epsilonMiddle, numE);
            logger.info() << "Mgr calculated for: [" << normalizedEnergies[bandIndices.front()] << ", ";
            logger << normalizedEnergies[bandIndices.back()] << "]." << std::endl;
            return bandIndices;
        }
    } else if (std::holds_alternative<CDFRange>(this->range)) {
        const auto &cdfRange = std::get<CDFRange>(this->range);
        double relativeIndexStart = cdfRange.cdfMiddle - cdfRange.cdfMargin / 2;
        double relativeIndexEnd = cdfRange.cdfMiddle + cdfRange.cdfMargin / 2;
        auto indexStart = static_cast<std::size_t>(eigensystem.size() * relativeIndexStart);
        auto indexEnd = static_cast<std::size_t>(eigensystem.size() * relativeIndexEnd);
        Assert(indexEnd < eigensystem.size());
        Assert(indexEnd > indexStart);

        std::vector<std::size_t> bandIndices(indexEnd - indexStart);
        std::iota(bandIndices.begin(), bandIndices.end(), indexStart);

        logger.info() << "Mgr calculated for: [" << normalizedEnergies[bandIndices.front()] << ", ";
        logger << normalizedEnergies[bandIndices.back()] << ")." << std::endl;

        return bandIndices;
    } else {
        throw std::runtime_error("Internal error");
    }
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

MeanGapRatio::EpsilonRange::EpsilonRange(double epsilonMiddle, Margin epsilonMargin)
        : epsilonMiddle{epsilonMiddle}, epsilonMargin{epsilonMargin}
{
    if (std::holds_alternative<WidthMargin>(epsilonMargin)) {
        double margin = std::get<WidthMargin>(epsilonMargin);
        Expects(margin > 0);
        Expects(epsilonMiddle - margin / 2 > 0 && epsilonMiddle + margin / 2 < 1);
    }
}

MeanGapRatio::VectorRange::VectorRange(FockBasis::Vector middleVector, Margin epsilonMargin)
        : middleVector{std::move(middleVector)}, epsilonMargin{epsilonMargin}
{
    if (std::holds_alternative<WidthMargin>(epsilonMargin)) {
        double margin = std::get<WidthMargin>(epsilonMargin);
        Expects(margin > 0);
    }
}

MeanGapRatio::CDFRange::CDFRange(double cdfMiddle, double cdfMargin) : cdfMiddle(cdfMiddle), cdfMargin(cdfMargin) {
    Expects(cdfMargin > 0);
    Expects(cdfMiddle - cdfMargin/2 > 0 && cdfMiddle + cdfMargin/2 < 1);
}