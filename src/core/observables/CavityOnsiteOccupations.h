//
// Created by Piotr Kubala on 06/12/2020.
//

#ifndef MBL_ED_CAVITYONSITEOCCUPATIONS_H
#define MBL_ED_CAVITYONSITEOCCUPATIONS_H

#include <memory>

#include "core/PrimaryObservable.h"
#include "core/FockBasis.h"
#include "core/terms/CavityLongInteraction.h"

/**
 * @brief Onsite mean occupations multiplied by a cosine like in cavity long interaction for all sites.
 * @details It is defined as
 * \f[ \langle \hat{n}_i \rangle \cos(2\pi\beta i + \phi_0). \f]
 */
class CavityOnsiteOccupations : public PrimaryObservable {
private:
    std::size_t numOfSites{};
    std::vector<arma::vec> diagonalObservables;
    std::vector<double> n_i{};
    std::shared_ptr<FockBasis> fockBasis;
    std::shared_ptr<CavityLongInteraction> cavityLongInteractions;

public:
    CavityOnsiteOccupations() = default;
    CavityOnsiteOccupations(std::shared_ptr<FockBasis> fockBasis,
                            std::shared_ptr<CavityLongInteraction> cavityLongInteractions);

    /**
     * @brief Returns the names of fields in format n_i_cos, where i = 1, ..., (number of sites)
     */
    [[nodiscard]] std::vector<std::string> getHeader() const override;

    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForState(const arma::cx_vec &state) override;
};


#endif //MBL_ED_CAVITYONSITEOCCUPATIONS_H
