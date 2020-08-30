//
// Created by Piotr Kubala on 10/02/2020.
//

#ifndef MBL_ED_CAVITYCONSTANTS_H
#define MBL_ED_CAVITYCONSTANTS_H

#include <vector>
#include <iosfwd>
#include <ostream>

/**
 * @brief Some constants from prof. Zakrzewski's notes, used in LookupCavityZ2 and LookupCavityYZ.
 * @details It is a set of different realisations for random @a phi0. Check CavityConstants::Realisation.
 */
class CavityConstants {
public:
    /**
     * @brief SiteEntry are 3 constants calculated from @a phi0 from CavityConstants::Realisation for a concrete site.
     */
    struct SiteEntry {
        /**
         * @brief It is just \f$ \cos(2\pi\beta i) \f$ for i-th site (like in CavityLongInteraction)
         */
        double cosine{};

        /**
         * @brief \f$ \int dx \cos(kx)w_i(x)^2 \f$.
         * @details It is the exact formula, which becomes SiteEntry::cosine if Wannier function \f$ w_i \f$ becomes
         * Dirac delta in a given site.
         */
        double wannier{};

        /**
         * @brief \f$ \int dx \cos(kx)w_i(x)w_{i+1}(x) \f$
         */
        double y{};

        friend bool operator==(const SiteEntry &lhs, const SiteEntry &rhs);
        friend bool operator!=(const SiteEntry &lhs, const SiteEntry &rhs);
        friend std::ostream &operator<<(std::ostream &os, const SiteEntry &entry);
    };

    /**
     * @brief Realisation is a concrete @a phi0 and a few site entries (see CavityConstants::SiteEntry) for each site
     * calculated from this @a phi0.
     */
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
    /**
     * @brief Adds another realisation to the set (presumably for another phi0).
     * @brief The number of sites is determined by the first invokation of this method and the following ones must
     * have the same number of CavityConstants::SiteEntry -ies.
     */
    void addRealisation(const Realisation &realisation);

    std::size_t getNumberOfSites() const;
    std::size_t size() const { return this->realisations.size(); }
    bool empty() const { return this->realisations.empty(); }
    const_iterator begin() const { return this->realisations.begin(); }
    const_iterator end() const { return this->realisations.end(); }
    const Realisation &operator[](std::size_t index) const { return this->realisations[index]; }
};


#endif //MBL_ED_CAVITYCONSTANTS_H
