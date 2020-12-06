//
// Created by Piotr Kubala on 06/12/2020.
//

#ifndef MBL_ED_CAVITYONSITEOCCUPATIONSSQUARED_H
#define MBL_ED_CAVITYONSITEOCCUPATIONSSQUARED_H

#include <memory>

#include <armadillo>

#include "core/PrimaryObservable.h"
#include "evolution/SymmetricMatrix.h"
#include "core/FockBasis.h"
#include "core/terms/CavityLongInteraction.h"

/**
 * @brief Expected values of all combinations of 2 occupation operators <n_j n_j>.
 */
class CavityOnsiteOccupationsSquared : public PrimaryObservable {
private:
    std::size_t numOfSites{};
    SymmetricMatrix<arma::vec> diagonalObservables;
    SymmetricMatrix<double> n_iN_j;
    SymmetricMatrix<std::string> headerStrings;
    std::shared_ptr<FockBasis> fockBasis;
    std::shared_ptr<CavityLongInteraction> cavityLongInteractions;

public:
    CavityOnsiteOccupationsSquared() = default;
    explicit CavityOnsiteOccupationsSquared(std::shared_ptr<FockBasis> fockBasis,
                                            std::shared_ptr<CavityLongInteraction> cavityLongInteractions);

    /**
     * @brief Returns the names of fields in format n_iN_j, where i, j = 1, ..., (number of sites).
     * @details The order is dictated by SymmetricMatrix class, namely n_iN_j correspond to
     * SymmetricMatrix::operator()(i, j) and the order is taken from SymmetricMatrix iterator.
     */
    [[nodiscard]] std::vector<std::string> getHeader() const override;

    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForState(const arma::cx_vec &state) override;

    /**
     * @brief More natural access to values that using getValues().
     */
    [[nodiscard]] virtual const SymmetricMatrix<double> &getOccupationsSquared() const;
};


#endif //MBL_ED_CAVITYONSITEOCCUPATIONSSQUARED_H
