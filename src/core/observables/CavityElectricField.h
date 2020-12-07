//
// Created by Piotr Kubala on 06/12/2020.
//

#ifndef MBL_ED_CAVITYELECTRICFIELD_H
#define MBL_ED_CAVITYELECTRICFIELD_H

#include "core/SecondaryObservable.h"

/**
 * @brief An observable representing a cavity output electric field, so a linear combination of occupation operators
 * with cosines appearing in cavity long interactions term.
 * @details It is given by
 *
 * \f[ \sum_{i=1}^K \cos(2\pi\beta i+\phi_0) \hat{n}_i. \f]
 *
 * For notation see CavityLongInteraction class.
 */
class CavityElectricField : public SecondaryObservable {
private:
    double electricField{};

public:
    [[nodiscard]] std::vector<std::string> getHeader() const override;
    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) override;
};


#endif //MBL_ED_CAVITYELECTRICFIELD_H
