//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_CAVITYLONGINTERACTION_H
#define MBL_ED_CAVITYLONGINTERACTION_H


#include "utils/Assertions.h"
#include "simulation/DiagonalTerm.h"

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
