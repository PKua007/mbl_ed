//
// Created by Piotr Kubala on 10/02/2020.
//

#include <iterator>
#include <functional>

#include "CavityConstantsReader.h"

namespace {
    template<typename Iterator, typename BinaryPredicate>
    void skipRepeating(Iterator &begin, const Iterator &end, BinaryPredicate predicate) {
        if (begin == end)
            return;

        using T = decltype(*begin);
        const T &first = *begin;

        begin++;
        while (begin != end) {
            if (!predicate(first, *begin))
                return;
            begin++;
        }
    }
}

/**
 * @brief Helper class for a single row from the file. Can convert itselt to CavityConstants::SiteEntry.
 */
struct CavityConstantsReader::Row {
    double phi0{};
    double cosine{};
    double wannier{};
    double y{};

    operator CavityConstants::SiteEntry() const { return {cosine, wannier, y}; }
    bool operator==(const Row &other) const { return this->phi0 == other.phi0; }

    friend std::istream &operator>>(std::istream &in, Row &row) {
        return in >> row.phi0 >> row.cosine >> row.wannier >> row.y;
    }
};

CavityConstants CavityConstantsReader::load(std::istream &in) {
    std::vector<Row> rows = loadRows(in);

    CavityConstants cavityConstants;
    if (rows.empty())
        return cavityConstants;

    std::size_t numberOfSites = countNumberOfSites(rows);
    if (rows.size() % numberOfSites != 0)
        throw std::runtime_error("Last realisation has too few elements");

    auto firstSiteEntryIt = rows.begin();
    while (firstSiteEntryIt != rows.end()) {
        auto nonEqualIt = std::adjacent_find(firstSiteEntryIt, firstSiteEntryIt + numberOfSites,
                                             std::not2(std::equal_to<Row>{}));
        if (nonEqualIt != firstSiteEntryIt + numberOfSites)
            throw std::runtime_error("Different number of sites in realisations");

        cavityConstants.addRealisation({firstSiteEntryIt->phi0, {firstSiteEntryIt, firstSiteEntryIt + numberOfSites}});
        firstSiteEntryIt += numberOfSites;
    }

    return cavityConstants;
}

std::size_t CavityConstantsReader::countNumberOfSites(const std::vector<Row> &rows) {
    auto it = rows.begin();
    skipRepeating(it, rows.end(), std::equal_to<Row>{});

    return it - rows.begin();
}

std::vector<CavityConstantsReader::Row> CavityConstantsReader::loadRows(std::istream &in) {
    std::vector<Row> rows;
    std::copy(std::istream_iterator<Row>(in), std::istream_iterator<Row>{}, std::back_inserter(rows));
    return rows;
}
