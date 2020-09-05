//
// Created by Piotr Kubala on 10/02/2020.
//

#include <stdexcept>
#include <tuple>
#include <istream>
#include <iterator>

#include "utils/Assertions.h"
#include "CavityConstants.h"

void CavityConstants::addRealisation(const CavityConstants::Realisation &realisation) {
    Expects(!realisation.siteEntries.empty());
    if (!this->empty())
        if (realisation.siteEntries.size() != this->getNumberOfSites())
            throw std::invalid_argument("Different number of site in the new entry");

    this->realisations.push_back(realisation);
}

std::size_t CavityConstants::getNumberOfSites() const {
    if (this->empty())
        return 0;
    return this->realisations.front().siteEntries.size();
}

bool operator==(const CavityConstants::SiteEntry &lhs, const CavityConstants::SiteEntry &rhs) {
    return std::tie(lhs.cosine, lhs.wannier, lhs.y) == std::tie(rhs.cosine, rhs.wannier, rhs.y);
}

bool operator!=(const CavityConstants::SiteEntry &lhs, const CavityConstants::SiteEntry &rhs) {
    return !(rhs == lhs);
}

std::ostream &operator<<(std::ostream &os, const CavityConstants::SiteEntry &entry) {
    return os << "{" << entry.cosine << ", " << entry.wannier << ", " << entry.y << "}";
}

bool operator==(const CavityConstants::Realisation &lhs, const CavityConstants::Realisation &rhs) {
    return std::tie(lhs.phi0, lhs.siteEntries) == std::tie(rhs.phi0, rhs.siteEntries);
}

bool operator!=(const CavityConstants::Realisation &lhs, const CavityConstants::Realisation &rhs) {
    return !(rhs == lhs);
}

std::ostream &operator<<(std::ostream &os, const CavityConstants::Realisation &realisation) {
    os << "{phi0: " << realisation.phi0 << ", siteEntries: {";
    std::copy(realisation.siteEntries.begin(), realisation.siteEntries.end(),
              std::ostream_iterator<CavityConstants::SiteEntry>(os, ", "));
    os << "}";
    return os;
}
