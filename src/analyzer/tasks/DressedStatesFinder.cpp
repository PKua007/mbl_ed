//
// Created by Piotr Kubala on 18/07/2020.
//

#include <algorithm>
#include <iterator>
#include <utility>

#include "DressedStatesFinder.h"

#include "utils/Assertions.h"

void DressedStatesFinder::analyze(const Eigensystem &eigensystem, Logger &logger) {
    Expects(eigensystem.hasEigenvectors());
    Expects(eigensystem.hasFockBasis());

    const auto &basis = eigensystem.getFockBasis();
    arma::mat normalizedEnergies = eigensystem.getNormalizedEigenenergies();
    auto indices = this->extractor.getBandIndices(eigensystem, logger);

    std::size_t numOfStatesFound{};
    for (std::size_t eigvecIdx : indices) {
        const auto &vector = eigensystem.getEigenstate(eigvecIdx);
        for (std::size_t coeffIndex{}; coeffIndex < vector.size(); coeffIndex++) {
            double coeffValue = vector[coeffIndex];
            if (std::abs(coeffValue) > this->coefficientThreshold) {
                std::ostringstream vectorStream;
                vectorStream << basis[coeffIndex];
                this->result.push_back({this->simulationIdx, vectorStream.str(), normalizedEnergies[eigvecIdx], coeffValue});
                numOfStatesFound++;
                break;
            }
        }
    }

    logger.info() << "Found " << numOfStatesFound << "/" << indices.size() << " dressed states. " << std::endl;
    this->simulationIdx++;
}

std::string DressedStatesFinder::getName() const {
    return "dressed";
}

void DressedStatesFinder::storeResult(std::ostream &out) const {
    std::copy(this->result.begin(), this->result.end(), std::ostream_iterator<Entry>(out, "\n"));
}

DressedStatesFinder::DressedStatesFinder(double coefficientThreshold, BandExtractor::Range range)
        : extractor(std::move(range), "Dressed states "), coefficientThreshold{coefficientThreshold}
{
    Expects(coefficientThreshold > M_SQRT1_2);
}

std::ostream &operator<<(std::ostream &out, const DressedStatesFinder::Entry &entry) {
    return out << entry.simulationIdx << " " << entry.vector << " " << entry.energy << " " << entry.coefficient;
}

void DressedStatesFinder::Entry::storeState(std::ostream &binaryOut) const {
    binaryOut.write(reinterpret_cast<const char*>(&this->simulationIdx), sizeof(this->simulationIdx));
    binaryOut.write(reinterpret_cast<const char*>(&this->energy), sizeof(this->energy));
    binaryOut.write(reinterpret_cast<const char*>(&this->coefficient), sizeof(this->coefficient));
    binaryOut.write(reinterpret_cast<const char*>(this->vector.c_str()), this->vector.size());
    binaryOut.write("\0", sizeof(char));
    Assert(binaryOut.good());
}

void DressedStatesFinder::Entry::restoreState(std::istream &binaryIn) {
    binaryIn.read(reinterpret_cast<char*>(&this->simulationIdx), sizeof(this->simulationIdx));
    binaryIn.read(reinterpret_cast<char*>(&this->energy), sizeof(this->energy));
    binaryIn.read(reinterpret_cast<char*>(&this->coefficient), sizeof(this->coefficient));
    std::getline(binaryIn, this->vector, '\0');
    Assert(binaryIn.good());
}

void DressedStatesFinder::storeState(std::ostream &binaryOut) const {
    binaryOut.write(reinterpret_cast<const char*>(&this->simulationIdx), sizeof(this->simulationIdx));
    std::size_t resultSize = this->result.size();
    binaryOut.write(reinterpret_cast<const char*>(&resultSize), sizeof(resultSize));
    Assert(binaryOut.good());
    for (const auto &entry : this->result)
        entry.storeState(binaryOut);
}

void DressedStatesFinder::joinRestoredState(std::istream &binaryIn) {
    std::size_t numOfSimulationsRestored{};
    binaryIn.read(reinterpret_cast<char*>(&numOfSimulationsRestored), sizeof(numOfSimulationsRestored));
    Assert(binaryIn.good());

    std::size_t resultSizeRestored{};
    binaryIn.read(reinterpret_cast<char*>(&resultSizeRestored), sizeof(resultSizeRestored));
    Assert(binaryIn.good());
    std::vector<Entry> entriesRestored(resultSizeRestored);
    for (auto &entry : entriesRestored) {
        entry.restoreState(binaryIn);
        entry.simulationIdx += this->simulationIdx;
    }
    this->result.reserve(this->result.size() + resultSizeRestored);
    this->result.insert(this->result.end(), entriesRestored.begin(), entriesRestored.end());

    this->simulationIdx += numOfSimulationsRestored;
}

void DressedStatesFinder::clear() {
    this->simulationIdx = 0;
    this->result.clear();
}
