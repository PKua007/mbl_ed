//
// Created by pkua on 01.11.2019.
//

#include <cmath>

#include "FockBase.h"
#include "Assertions.h"

void FockBase::add(FockBase::Vector vector) {
    if (this->size() != 0)
        Expects(this->theBase.front().size() == vector.size());

    std::size_t currIndex = this->size();
    this->indexMap.emplace(this->computeHash(vector), currIndex);
    this->theBase.push_back(std::move(vector));
}

std::size_t FockBase::size() const {
    return this->theBase.size();
}

FockBase::Vector &FockBase::operator[](std::size_t i) {
    return this->theBase[i];
}

const FockBase::Vector &FockBase::operator[](std::size_t i) const {
    return this->theBase[i];
}

std::optional<std::size_t> FockBase::findIndex(const FockBase::Vector &vector) const {
    auto it = this->indexMap.find(this->computeHash(vector));
    if (it == this->indexMap.end())
        return std::nullopt;
    else
        return it->second;
}

FockBase::iterator FockBase::begin() {
    return this->theBase.begin();
}

FockBase::iterator FockBase::end() {
    return this->theBase.end();
}

FockBase::const_iterator FockBase::begin() const {
    return this->theBase.begin();
}

FockBase::const_iterator FockBase::end() const {
    return this->theBase.end();
}

double FockBase::computeHash(const FockBase::Vector &vector) const {
    double hash{};
    for (std::size_t i = 0; i < vector.size(); i++)
        hash += std::sqrt(100*i + 3) * vector[i];
    return hash;
}
