//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_HUBBARDHOP_H
#define MBL_ED_HUBBARDHOP_H


#include "utils/Assertions.h"
#include "simulation/HoppingTerm.h"

/**
 * @brief A standard one-site hop from Hubbard Hamiltonian.
 * @details It is defined as
 * \f[ -J \sum_{i=1}^K ( \hat{b}_i^\dagger \hat{b}_{i+1} + \text{c. c.} ) \f]
 *
 * where \f$ i \f$ is the number of site, \f$ K \f$ is the total number of sites, \f$ \hat{b}_i \f$ is annihilation
 * operator and the constant \f$ J \f$ is passed in the constructor.
 */
class HubbardHop : public HoppingTerm {
private:
    double J{};

public:
    explicit HubbardHop(double J);

    /**
     * @brief Returns constant -J for a hop between neighbouring sites (according to PBC or OBC).
     */
    double calculate(const FockBase::Vector &from, const FockBase::Vector &to, std::size_t fromSite,
                     std::size_t toSite, const HamiltonianGenerator &generator) override;
};

#endif //MBL_ED_HUBBARDHOP_H
