//
// Created by Piotr Kubala on 06/03/2021.
//

#include <cmath>

#include "ParticipationEntropy.h"
#include "simulation/RestorableHelper.h"

ParticipationEntropy::ParticipationEntropy(double q, BandExtractor::Range range)
        : extractor(std::move(range), "Participation entropy"), q{q}
{
    Expects(q > 0);
}

void ParticipationEntropy::analyze(const Eigensystem &eigensystem, [[maybe_unused]] Logger &logger) {
    Expects(eigensystem.hasEigenvectors());
    auto normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    auto bandIndices = this->extractor.getBandIndices(eigensystem, logger);

    double singleEntropy{};
    std::size_t numEntries{};
    for (auto i : bandIndices) {
        auto state = eigensystem.getEigenstate(i);
        double entropy{};
        if (this->q == 1) {     // Shannon's entropy limit
            entropy = std::accumulate(state.begin(), state.end(), 0., [](double sum, double d) {
                if (d == 0)
                    return sum;
                else
                    return sum - d*d * std::log(d*d);
            });
        } else {    // General case
            entropy = std::accumulate(state.begin(), state.end(), 0., [this](double sum, double d) {
                return sum + std::pow(d*d, this->q);
            });
            entropy = 1/(1 - this->q) * std::log(entropy);
        }
        singleEntropy += entropy;
        numEntries++;
    }
    if (numEntries > 0)
        this->entropies.push_back(singleEntropy / numEntries);
}

std::string ParticipationEntropy::getName() const {
    return "pe";
}

std::vector<std::string> ParticipationEntropy::getResultHeader() const {
    return {"participationEntropy", "participationEntropyError"};
}

std::vector<std::string> ParticipationEntropy::getResultFields() const {
    Quantity result = this->calculateMean();
    result.separator = Quantity::Separator::SPACE;
    std::stringstream resultStream;
    resultStream << result;
    std::string quantityValue, quantityError;
    resultStream >> quantityValue >> quantityError;
    return {quantityValue, quantityError};
}

void ParticipationEntropy::storeState(std::ostream &binaryOut) const {
    RestorableHelper::storeStateForVector(this->entropies, binaryOut);
}

void ParticipationEntropy::joinRestoredState(std::istream &binaryIn) {
    RestorableHelper::joinRestoredStateForVector(this->entropies, binaryIn);
}

void ParticipationEntropy::clear() {
    this->entropies.clear();
}

Quantity ParticipationEntropy::calculateMean() const {
    Quantity result;
    result.calculateFromSamples(this->entropies);
    return result;
}
