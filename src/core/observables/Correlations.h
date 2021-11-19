//
// Created by Piotr Kubala on 21/09/2020.
//

#ifndef MBL_ED_CORRELATIONS_H
#define MBL_ED_CORRELATIONS_H

#include "core/SecondaryObservable.h"

/**
 * @brief Site-averaged correlations for all sites apart from @a marginSize sites on the border.
 * @details Correlations are defined as G(i, i+d) = \<n_i n_{i+d}\> - \<n_i\> \<n_{i+d}\>. The class outputs
 * correlations for all d-s which are plausible - namely at least one valule will go to site averaging. As an example,
 * for 6 sites and margin 1, we have: G_1 as an average of G(2, 3), G(3, 4), G(4, 5); G_2 - G(2, 4), G(3, 5);
 * G_3 - G(2, 5).
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

    /**
     * @brief Returns the names of fields with correlations, names Gmx_d, where 'm' is margin size and 'd' is described
     * the class description.
     */
    [[nodiscard]] std::vector<std::string> getHeader() const override;

    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) override;
};

#endif //MBL_ED_CORRELATIONS_H
