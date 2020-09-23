//
// Created by Piotr Kubala on 21/09/2020.
//

#ifndef MBL_ED_CORRELATIONS_H
#define MBL_ED_CORRELATIONS_H

#include "evolution/SecondaryObservable.h"

/**
 * @brief Site-averaged correlations for all sites apart from @a marginSize sites on the border.
 * @details Correlations are defined as G(i, i+d) = <n_i n_{i+d}> - <n_i> <n_{i+d}>. The struct contains
 * correlations for a specific @a d - distance, averaged over all sites i without the before mentioned marginal
 * sites. So, for example for 5 sites, marginSize = 1, d = 1, G_d = (G(1, 2) + G(2, 3)) / 2. The values are
 * calculated from OccupationEvolution::Occupations observables and can be added multiple times - they are then
 * averaged.
 */
class Correlations : public SecondaryObservable {
private:
    std::size_t numOfSites{};
    std::size_t marginSize{};
    std::size_t borderlessNumOfSites{};
    std::vector<double> G_d{};

public:
    Correlations() = default;
    Correlations(std::size_t numOfSites, std::size_t marginSize);

    [[nodiscard]] std::vector<std::string> getHeader() const override;
    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) override;
};

#endif //MBL_ED_CORRELATIONS_H
