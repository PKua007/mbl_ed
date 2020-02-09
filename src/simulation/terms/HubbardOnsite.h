//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_HUBBARDONSITE_H
#define MBL_ED_HUBBARDONSITE_H


#include "utils/Assertions.h"
#include "simulation/DiagonalTerm.h"

class HubbardOnsite : public DiagonalTerm {
private:
    double U{};

public:
    explicit HubbardOnsite(double U);

    double calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) override;
};


#endif //MBL_ED_HUBBARDONSITE_H
