//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_CAVITYLONGINTERACTION_H
#define MBL_ED_CAVITYLONGINTERACTION_H


#include "utils/Assertions.h"
#include "core/DiagonalTerm.h"

/**
 * @brief A cavity-mediated long range interactions.
 * @details It is defined as
 * \f[ -\frac{U_1}{K} \left( \sum_{i=1}^K \cos(2\pi\beta i+\phi_0) \hat{n}_i \right)^2, \f]
 *
 * where \f$ i \f$ is the number of site, \f$ K \f$ is the total number of sites, \f$ \hat{b}_i \f$ is annihilation
 * operator, \f$ \hat{n}_i = \hat{b}_i^\dagger\hat{b}_i \f$, and the rest are the parameters passed in the constructor.
 */
class CavityLongInteraction : public DiagonalTerm {
private:
    double U1{};
    double beta{};
    double phi0{};

public:
    CavityLongInteraction(double U1, double beta, double phi0);

    double calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) override;

    void setPhi0(double phi0);
};


#endif //MBL_ED_CAVITYLONGINTERACTION_H
