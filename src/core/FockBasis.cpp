//
// Created by pkua on 01.11.2019.
//

#include <cmath>
#include <numeric>

#include "FockBasis.h"
#include "utils/Assertions.h"

void FockBasis::add(FockBasis::Vector vector) {
    if (this->size() != 0)
        Expects(this->theBasis.front().size() == vector.size());

    std::size_t currIndex = this->size();
    this->indexMap.emplace(this->computeHash(vector), currIndex);
    this->theBasis.push_back(std::move(vector));
}

std::size_t FockBasis::size() const {
    return this->theBasis.size();
}

FockBasis::Vector &FockBasis::operator[](std::size_t i) {
    return this->theBasis[i];
}

const FockBasis::Vector &FockBasis::operator[](std::size_t i) const {
    return this->theBasis[i];
}

std::optional<std::size_t> FockBasis::findIndex(const FockBasis::Vector &vector) const {
    if (vector.size() != this->getNumberOfSites())
        return std::nullopt;

    auto it = this->indexMap.find(this->computeHash(vector));
    if (it == this->indexMap.end())
        return std::nullopt;
    else
        return it->second;
}

FockBasis::iterator FockBasis::begin() {
    return this->theBasis.begin();
}

FockBasis::iterator FockBasis::end() {
    return this->theBasis.end();
}

FockBasis::const_iterator FockBasis::begin() const {
    return this->theBasis.begin();
}

FockBasis::const_iterator FockBasis::end() const {
    return this->theBasis.end();
}

double FockBasis::computeHash(const FockBasis::Vector &vector) const {
    // https://arxiv.org/pdf/1102.4006.pdf says, that sqrt(100*i + 3) is linearly independent over radicals, so
    // each distinct vector will end up with a unique hash
    double hash{};
    for (std::size_t i = 0; i < vector.size(); i++)
        hash += std::sqrt(100*i + 3) * vector[i];
    return hash;
}

std::size_t FockBasis::getNumberOfSites() const {
    Expects(this->size() > 0);
    return this->theBasis.front().size();
}

std::size_t FockBasis::getNumberOfParticles() const {
    Expects(this->size() > 0);
    return std::accumulate(this->theBasis.front().begin(), this->theBasis.front().end(), 0);
}
