//
// Created by Piotr Kubala on 21/09/2020.
//

#ifndef MBL_ED_ONSITEOCCUPATIONS_H
#define MBL_ED_ONSITEOCCUPATIONS_H

#include <armadillo>

#include "evolution/PrimaryObservable.h"
#include "core/FockBase.h"

/**
 * @brief Onsite mean occupations for a site @a i defined as <n_i>. The values are calculated from
 * OccupationEvolution::Occupations observables and can be added multiple times - they are then averaged.
 */
class OnsiteOccupations : public PrimaryObservable {
private:
    std::size_t numOfSites{};
    std::vector<arma::vec> diagonalObservables;
    std::vector<double> n_i{};
    std::shared_ptr<FockBase> fockBase;

public:
    OnsiteOccupations() = default;
    OnsiteOccupations(std::size_t numOfSites, std::shared_ptr<FockBase> fockBase);

    [[nodiscard]] std::vector<std::string> getHeader() const override;
    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForState(const arma::cx_vec &state) override;
};


#endif //MBL_ED_ONSITEOCCUPATIONS_H
