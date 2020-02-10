//
// Created by Piotr Kubala on 10/02/2020.
//

#ifndef MBL_ED_CAVITYCONSTANTS_H
#define MBL_ED_CAVITYCONSTANTS_H

#include <vector>
#include <iosfwd>
#include <ostream>

class CavityConstants {
public:
    struct SiteEntry {
        double cosine{};
        double wannier{};
        double y{};

        friend bool operator==(const SiteEntry &lhs, const SiteEntry &rhs);
        friend bool operator!=(const SiteEntry &lhs, const SiteEntry &rhs);
        friend std::ostream &operator<<(std::ostream &os, const SiteEntry &entry);
    };

    struct Realisation {
        double phi0{};
        std::vector<SiteEntry> siteEntries;

        friend bool operator==(const Realisation &lhs, const Realisation &rhs);
        friend bool operator!=(const Realisation &lhs, const Realisation &rhs);

        friend std::ostream &operator<<(std::ostream &os, const Realisation &realisation);
    };

    using iterator = std::vector<Realisation>::iterator;
    using const_iterator = std::vector<Realisation>::const_iterator;

private:
    std::vector<Realisation> realisations;

public:
    void addRealisation(const Realisation &realisation);

    std::size_t getNumberOfSites() const;
    std::size_t size() const { return this->realisations.size(); }
    bool empty() const { return this->realisations.empty(); }
    const_iterator begin() const { return this->realisations.begin(); }
    const_iterator end() const { return this->realisations.end(); }
    const Realisation &operator[](std::size_t index) const { return this->realisations[index]; }
};


#endif //MBL_ED_CAVITYCONSTANTS_H
