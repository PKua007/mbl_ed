//
// Created by Piotr Kubala on 21/09/2020.
//

#include <sstream>
#include <iterator>

#include "TimeEvolutionParameters.h"
#include "utils/Assertions.h"


void TimeEvolutionParameters::setVectorsToEvolveFromTags(const std::vector <std::string> &strings) {
    Expects(this->numberOfSites > 0);

    this->vectorsToEvolve.clear();
    for (const auto &string : strings) {
        try {
            // Try a tag representation
            this->vectorsToEvolve.emplace_back(FockBasis::Vector(this->numberOfSites, string));
        } catch (FockVectorParseException &) {
            // If tag failed, try occupation representation
            this->vectorsToEvolve.emplace_back(FockBasis::Vector(string));
        }
    }
}

std::size_t TimeEvolutionParameters::countStoredObservableValues() const {
    std::size_t numOfValues{};
    for (const auto &storedObservable : this->storedObservables)
        numOfValues += storedObservable->getHeader().size();
    return numOfValues;
}

std::string TimeEvolutionParameters::generateStoredObservablesHeader() const {
    std::vector<std::string> headerStrings;
    headerStrings.reserve(this->countStoredObservableValues());
    for (const auto &storedObservable : this->storedObservables) {
        auto observableHeaderStrings = storedObservable->getHeader();
        headerStrings.insert(headerStrings.end(), observableHeaderStrings.begin(), observableHeaderStrings.end());
    }

    std::ostringstream out;
    std::copy(headerStrings.begin(), headerStrings.end(), std::ostream_iterator<std::string>(out, " "));
    return out.str();
}
