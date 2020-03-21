//
// Created by Piotr Kubala on 16/03/2020.
//

#ifndef MBL_ED_QUASIPERIODICDISORDER_H
#define MBL_ED_QUASIPERIODICDISORDER_H


#include "simulation/DiagonalTerm.h"

/**
 * @brief The class representing quasiperiodic potential in each site.
 * @details It is defined as
 * \f[ \sum_{i=1}^K W \cos(2\pi \beta i + \phi_0) \hat{n}_i \f]
 *
 * where \f$ i \f$ is the number of site, \f$ K \f$ is the total number of sites, \f$ \hat{b}_i \f$ is annihilation
 * operator, \f$ \hat{n}_i = \hat{b}_i^\dagger\hat{b}_i \f$. The rest of parameters are in the constructor (also
 * \f$ \phi_0 \f$ can be changed by setPhi0().
 * @tparam DisorderGenerator The concrete disorded generator, whose operator() is used to sample the disorder.
 */
class QuasiperiodicDisorder : public DiagonalTerm {double U1{};
    double W{};
    double beta{};
    double phi0{};

public:
    QuasiperiodicDisorder(double W, double beta, double phi0);

    double calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) override;

    void setPhi0(double phi0);
};


#endif //MBL_ED_QUASIPERIODICDISORDER_H
