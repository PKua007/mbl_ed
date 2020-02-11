//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_HUBBARDONSITE_H
#define MBL_ED_HUBBARDONSITE_H


#include "utils/Assertions.h"
#include "simulation/DiagonalTerm.h"

/**
 * @brief A stantard onsite energy term from Hubbard hamiltonian.
 * @details It is defined as
 * \f[ \frac{U}{2} \sum_{i=1}^K \hat{n}_i (\hat{n}_i - 1) \f]
 *
 * where \f$ i \f$ is the number of site, \f$ K \f$ is the total number of sites, \f$ \hat{b}_i \f$ is annihilation
 * operator, \f$ \hat{n}_i = \hat{b}_i^\dagger\hat{b}_i \f$ and \f$ U \f$ is a parameter in the constructor.
 */
class HubbardOnsite : public DiagonalTerm {
private:
    double U{};

public:
    explicit HubbardOnsite(double U);

    double calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) override;
};


#endif //MBL_ED_HUBBARDONSITE_H
