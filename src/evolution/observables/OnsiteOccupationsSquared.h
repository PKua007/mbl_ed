//
// Created by Piotr Kubala on 21/09/2020.
//

#ifndef MBL_ED_ONSITEOCCUPATIONSSQUARED_H
#define MBL_ED_ONSITEOCCUPATIONSSQUARED_H

#include <armadillo>

#include "evolution/PrimaryObservable.h"
#include "evolution/SymmetricMatrix.h"
#include "core/FockBase.h"

/**
 * @brief Onsite mean occupations for a site @a i defined as <n_i>. The values are calculated from
 * OccupationEvolution::Occupations observables and can be added multiple times - they are then averaged.
 */
class OnsiteOccupationsSquared : public PrimaryObservable {
private:
    std::size_t numOfSites{};
    SymmetricMatrix<arma::vec> diagonalObservables;
    SymmetricMatrix<double> n_iN_j;
    SymmetricMatrix<std::string> headerStrings;
    std::shared_ptr<FockBase> fockBase;

public:
    OnsiteOccupationsSquared() = default;
    explicit OnsiteOccupationsSquared(std::shared_ptr<FockBase> fockBase);

    [[nodiscard]] std::vector<std::string> getHeader() const override;
    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForState(const arma::cx_vec &state) override;

    virtual const SymmetricMatrix<double> &getOccupationsSquared() const;
};

#endif //MBL_ED_ONSITEOCCUPATIONSSQUARED_H
