//
// Created by Piotr Kubala on 21/09/2020.
//

#ifndef MBL_ED_ONSITEFLUCTUATIONS_H
#define MBL_ED_ONSITEFLUCTUATIONS_H

#include "core/SecondaryObservable.h"

/**
 * @brief Onsite fluctuations for all sites defined as rho(i) = \<n_i^2\> - \<n_i\>^2.
 */
class OnsiteFluctuations : public SecondaryObservable {
private:
    std::size_t numOfSites{};
    std::vector<double> rho_i{};

public:
    OnsiteFluctuations() = default;
    explicit OnsiteFluctuations(std::size_t numOfSites);

    /**
     * @brief Returns the names of fields in format rho_i, where i = 1, ..., (number of sites)
     */
    [[nodiscard]] std::vector<std::string> getHeader() const override;

    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) override;
};


#endif //MBL_ED_ONSITEFLUCTUATIONS_H
