//
// Created by Piotr Kubala on 06/12/2020.
//

#ifndef MBL_ED_CAVITYLIGHTINTENSITY_H
#define MBL_ED_CAVITYLIGHTINTENSITY_H


#include "core/SecondaryObservable.h"

/**
 * @brief An observable representing a cavity output light intensity, so basically the same term as in cavity mediated
 * long range interactions
 * @details It is given by
 *
 * \f[ \left( \sum_{i=1}^K \cos(2\pi\beta i+\phi_0) \hat{n}_i \right)^2. \f]
 *
 * For notation see CavityLongInteraction class.
 */
class CavityLightIntensity : public SecondaryObservable {
private:
    double lightIntensity{};

public:
    [[nodiscard]] std::vector<std::string> getHeader() const override;
    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) override;
};


#endif //MBL_ED_CAVITYLIGHTINTENSITY_H
