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

    struct Row {
        double phi0{};
        double cosine{};
        double wannier{};
        double y{};

        operator CavityConstants::SiteEntry() const { return {cosine, wannier, y}; }
        bool operator==(const Row &other) const { return this->phi0 == other.phi0; }
    };

    std::istream &operator>>(std::istream &in, Row &row) {
        return in >> row.phi0 >> row.cosine >> row.wannier >> row.y;
    }
}

CavityConstants CavityConstantsReader::load(std::istream &in) {
    std::vector<Row> entries;
    std::copy(std::istream_iterator<Row>(in), std::istream_iterator<Row>{}, std::back_inserter(entries));

    CavityConstants cavityConstants;
    if (entries.empty())
        return cavityConstants;

    auto it = entries.begin();
    skipRepeating(it, entries.end(), std::equal_to<Row>{});

    std::size_t numberOfSites = it - entries.begin();
    if (entries.size() % numberOfSites != 0)
        throw std::runtime_error("Last realisation has too few elements");

    auto it2 = entries.begin();
    while (it2 != entries.end()) {
        auto nonEqualIt = std::adjacent_find(it2, it2 + numberOfSites, std::not2(std::equal_to<Row>{}));
        if (nonEqualIt != it2 + numberOfSites)
            throw std::runtime_error("Different number of sites in realisations");

        cavityConstants.addRealisation({it2->phi0, {it2, it2 + numberOfSites}});
        it2 += numberOfSites;
    }

    return cavityConstants;
}
