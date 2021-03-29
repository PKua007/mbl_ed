//
// Created by Piotr Kubala on 28/03/2021.
//

#include "BandExtractor.h"
#include "utils/Assertions.h"

BandExtractor::EpsilonRange::EpsilonRange(double epsilonMiddle, Margin epsilonMargin)
        : epsilonMiddle{epsilonMiddle}, epsilonMargin{epsilonMargin}
{
    if (std::holds_alternative<WidthMargin>(epsilonMargin)) {
        double margin = std::get<WidthMargin>(epsilonMargin);
        Expects(margin > 0);
        Expects(epsilonMiddle - margin / 2 > 0 && epsilonMiddle + margin / 2 < 1);
    }
}

BandExtractor::VectorRange::VectorRange(FockBasis::Vector middleVector, Margin epsilonMargin)
        : middleVector{std::move(middleVector)}, epsilonMargin{epsilonMargin}
{
    if (std::holds_alternative<WidthMargin>(epsilonMargin)) {
        double margin = std::get<WidthMargin>(epsilonMargin);
        Expects(margin > 0);
    }
}

BandExtractor::CDFRange::CDFRange(double cdfMiddle, double cdfMargin)
        : cdfMiddle(cdfMiddle), cdfMargin(cdfMargin)
{
    Expects(cdfMargin > 0);
    Expects(cdfMiddle - cdfMargin/2 > 0 && cdfMiddle + cdfMargin/2 < 1);
}

std::vector<std::size_t> BandExtractor::getBandIndices(const Eigensystem &eigensystem, Logger &logger) const {
    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();

    if (std::holds_alternative<VectorRange>(this->range))
        return this->getIndicesForVectorRange(eigensystem, normalizedEnergies, std::get<VectorRange>(range), logger);
    else if (std::holds_alternative<EpsilonRange>(this->range))
        return this->getIndicesForEpsilonRange(eigensystem, normalizedEnergies, std::get<EpsilonRange>(range), logger);
    else if (std::holds_alternative<CDFRange>(this->range))
        return this->getIndicesForCDFRange(eigensystem, normalizedEnergies, std::get<CDFRange>(range), logger);
    else
        throw std::runtime_error("Internal error");
}

std::vector<std::size_t> BandExtractor::getIndicesForCDFRange(const Eigensystem &eigensystem,
                                                              const arma::vec &normalizedEnergies,
                                                              const BandExtractor::CDFRange &cdfRange,
                                                              Logger &logger) const
{
    double relativeIndexStart = cdfRange.cdfMiddle - cdfRange.cdfMargin / 2;
    double relativeIndexEnd = cdfRange.cdfMiddle + cdfRange.cdfMargin / 2;
    auto indexStart = static_cast<std::size_t>(eigensystem.size() * relativeIndexStart);
    auto indexEnd = static_cast<std::size_t>(eigensystem.size() * relativeIndexEnd);
    Assert(indexEnd < eigensystem.size());
    Assert(indexEnd > indexStart);

    std::vector<std::size_t> bandIndices(indexEnd - indexStart);
    std::iota(bandIndices.begin(), bandIndices.end(), indexStart);

    logger.info() << this->taskName << " calculated for: [" << normalizedEnergies[bandIndices.front()] << ", ";
    logger << normalizedEnergies[bandIndices.back()] << ")." << std::endl;

    return bandIndices;
}

std::vector<std::size_t> BandExtractor::getIndicesForEpsilonRange(const Eigensystem &eigensystem,
                                                                  const arma::vec &normalizedEnergies,
                                                                  const BandExtractor::EpsilonRange &epsilonRange,
                                                                  Logger &logger) const
{
    if (std::holds_alternative<WidthMargin>(epsilonRange.epsilonMargin)) {
        double margin = std::get<WidthMargin>(epsilonRange.epsilonMargin);
        return eigensystem.getIndicesOfNormalizedEnergiesInBand(epsilonRange.epsilonMiddle, margin);
    } else { // holds alternative NumberOfEnergiesMargin
        std::size_t numEnergies = std::get<NumberOfEnergiesMargin>(epsilonRange.epsilonMargin);
        auto bandIndices = eigensystem.getIndicesOfNumberOfNormalizedEnergies(epsilonRange.epsilonMiddle, numEnergies);
        logger.info() << this->taskName << " calculated for: [" << normalizedEnergies[bandIndices.front()] << ", ";
        logger << normalizedEnergies[bandIndices.back()] << "]." << std::endl;
        return bandIndices;
    }
}

std::vector<std::size_t> BandExtractor::getIndicesForVectorRange(const Eigensystem &eigensystem,
                                                                 const arma::vec &normalizedEnergies,
                                                                 const BandExtractor::VectorRange &vectorRange,
                                                                 Logger &logger) const
{
    double relativeMiddleEnergy = calculateEnergyOfFockState(vectorRange.middleVector, eigensystem);

    if (std::holds_alternative<WidthMargin>(vectorRange.epsilonMargin)) {
        double relativeMargin = std::get<WidthMargin>(vectorRange.epsilonMargin);
        if (relativeMiddleEnergy - relativeMargin / 2 <= 0 || relativeMiddleEnergy + relativeMargin / 2 >= 1) {
            throw std::runtime_error(
                "Margin " + std::to_string(relativeMargin) + " is to big around the given vector of the energy "
                + std::to_string(relativeMiddleEnergy)
            );
        }

        logger.info() << this->taskName << " calculated around: " << relativeMiddleEnergy << ". " << std::endl;
        return eigensystem.getIndicesOfNormalizedEnergiesInBand(relativeMiddleEnergy, relativeMargin);
    } else { // holds alternative NumberOfEnergiesMargin
        std::size_t numEnergies = std::get<NumberOfEnergiesMargin>(vectorRange.epsilonMargin);
        auto bandIndices = eigensystem.getIndicesOfNumberOfNormalizedEnergies(relativeMiddleEnergy, numEnergies);
        logger.info() << this->taskName << " calculated around: " << relativeMiddleEnergy << ", for: [";
        logger << normalizedEnergies[bandIndices.front()] << ", " << normalizedEnergies[bandIndices.back()] << "].";
        logger << std::endl;
        return bandIndices;
    }
}

double BandExtractor::calculateEnergyOfFockState(const FockBasis::Vector &state,
                                                 const Eigensystem &eigensystem)
{
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
