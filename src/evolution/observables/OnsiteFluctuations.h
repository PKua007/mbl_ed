//
// Created by Piotr Kubala on 21/09/2020.
//

#ifndef MBL_ED_ONSITEFLUCTUATIONS_H
#define MBL_ED_ONSITEFLUCTUATIONS_H

#include "evolution/SecondaryObservable.h"

/**
 * @brief Onsite fluctuations for a site @a i defined as rho(i) = <n_i^2> - <n_i>^2. The values are calculated from
 * OccupationEvolution::Occupations observables and can be added multiple times - they are then averaged.
 */
class OnsiteFluctuations : public SecondaryObservable {
private:
    std::size_t numOfSites{};
    std::vector<double> rho_i{};

public:
    OnsiteFluctuations() = default;
    explicit OnsiteFluctuations(std::size_t numOfSites);

    [[nodiscard]] std::vector<std::string> getHeader() const override;
    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) override;
};


#endif //MBL_ED_ONSITEFLUCTUATIONS_H
