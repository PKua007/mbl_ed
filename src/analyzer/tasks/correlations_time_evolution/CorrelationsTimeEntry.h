//
// Created by Piotr Kubala on 21/02/2020.
//

#ifndef MBL_ED_CORRELATIONSTIMEENTRY_H
#define MBL_ED_CORRELATIONSTIMEENTRY_H

#include <vector>

#include "OccupationEvolution.h"

/**
 * @brief Class collecting OccupationEvolution::Occupations for subsequent eigensystems and calculating disorder
 * averages of correlations as well as onsite occupations and fluctuations calculated from before mentioned occupations
 * for a specific moment of time.
 * @details The class contains:
 * <ul>
 *     <li>site-averaged correlations, for all distanced from 1 to <em> number of sites - 1 </em>
 *     (see CorrelationsTimeEntry::Correlations in the code)
 *     <li>site-averaget correlations, but taking into account only sites after removing a few on the borders
 *     (@a marginSize parameters from the constructor)
 *     <li>onsite fluctuations (see CorrelationsTimeEntry::OnsiteFluctuations in the code)
 *     <li>onsite occupaitons (see CorrelationsTimeEntry::OnsiteOccupations in the code)
 * </ul>
 */
class CorrelationsTimeEntry {
private:
    /**
     * @brief Site-averaged correlations for all sites apart from @a marginSize sites on the border.
     * @details Correlations are defined as G(i, i+d) = <n_i n_{i+d}> - <n_i> <n_{i+d}>. The struct contains
     * correlations for a specific @a d - distance, averaged over all sites i without the before mentioned marginal
     * sites. So, for example for 5 sites, marginSize = 1, d = 1, G_d = (G(1, 2) + G(2, 3)) / 2. The values are
     * calculated from OccupationEvolution::Occupations observables and can be added multiple times - they are then
     * averaged.
     */
    class Correlations {
    private:
        std::size_t d{};
        std::size_t marginSize{};
        double G_d{};
        std::size_t numberOfMeanEntries{};

    public:
        Correlations() = default;
        Correlations(size_t d, size_t marginSize) : d{d}, marginSize{marginSize} { }

        void addObservables(const OccupationEvolution::Occupations &occupations);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string toString() const;
    };

    /**
     * @brief Onsite fluctuations for a site @a i defined as rho(i) = <n_i^2> - <n_i>^2. The values are calculated from
     * OccupationEvolution::Occupations observables and can be added multiple times - they are then averaged.
     */
    class OnsiteFluctuations {
    private:
        std::size_t i{};
        double rho_i{};
        std::size_t numberOfMeanEntries{};

    public:
        OnsiteFluctuations() = default;
        explicit OnsiteFluctuations(std::size_t i) : i{i} { }

        void addObservables(const OccupationEvolution::Occupations &occupations);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string toString() const;
    };

    /**
     * @brief Onsite mean occupations for a site @a i defined as <n_i>. The values are calculated from
     * OccupationEvolution::Occupations observables and can be added multiple times - they are then averaged.
     */
    class OnsiteOccupations {
    private:
        std::size_t i{};
        double n_i{};
        std::size_t numberOfMeanEntries{};

    public:
        OnsiteOccupations() = default;
        explicit OnsiteOccupations(std::size_t i) : i{i} { }

        void addObservables(const OccupationEvolution::Occupations &occupations);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string toString() const;
    };

    double t{};
    std::vector<Correlations> correlations{};
    std::vector<Correlations> borderlessCorrelations{};
    std::vector<OnsiteFluctuations> onsiteFluctuations{};
    std::vector<OnsiteOccupations> onsiteOccupations{};

    std::size_t numberOfSites{};
    std::size_t numberOfMeanEntries{};

    void populateOnsiteFluctuations(std::vector<OnsiteFluctuations> &onsiteFluctuationsVector) const;
    void populateOnsiteOccupations(std::vector<OnsiteOccupations> &onsiteOccupationsVector) const;
    void populateCorrelations(std::vector<Correlations> &correlationsVector, std::size_t marginSize_) const;

public:
    CorrelationsTimeEntry() = default;
    CorrelationsTimeEntry(double t, std::size_t marginSize, std::size_t numberOfSites);

    [[nodiscard]] std::size_t getNumberOfSites() const;

    /**
     * @brief Adds next entry of observables
     */
    void addObservables(const OccupationEvolution::Occupations &occupations);

    /**
     * @brief Constructs the header for all stuff - time, correlation, correlations without margin and fluctuations.
     */
    [[nodiscard]] std::string getHeader() const;

    /**
     * @brief Constructs the row of values of all stuff - time, correlation, correlations without margin and
     * fluctuations.
     */
    [[nodiscard]] std::string toString() const;
};


#endif //MBL_ED_CORRELATIONSTIMEENTRY_H
