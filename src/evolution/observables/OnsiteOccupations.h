//
// Created by Piotr Kubala on 21/09/2020.
//

#ifndef MBL_ED_ONSITEOCCUPATIONS_H
#define MBL_ED_ONSITEOCCUPATIONS_H

#include <memory>

#include <armadillo>

#include "evolution/PrimaryObservable.h"
#include "core/FockBasis.h"

/**
 * @brief Onsite mean occupations for all sites defined as <n_i>.
 */
class OnsiteOccupations : public PrimaryObservable {
private:
    std::size_t numOfSites{};
    std::vector<arma::vec> diagonalObservables;
    std::vector<double> n_i{};
    std::shared_ptr<FockBasis> fockBasis;

public:
    OnsiteOccupations() = default;
    OnsiteOccupations(std::shared_ptr<FockBasis> fockBasis);

    /**
     * @brief Returns the names of fields in format n_i, where i = 1, ..., (number of sites)
     */
    [[nodiscard]] std::vector<std::string> getHeader() const override;

    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForState(const arma::cx_vec &state) override;
};


#endif //MBL_ED_ONSITEOCCUPATIONS_H
