//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_DIAGONALTERM_H
#define MBL_ED_DIAGONALTERM_H

#include "FockBase.h"

class HamiltonianGenerator;

class DiagonalTerm {
public:
    virtual ~DiagonalTerm() = default;

    virtual double calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) = 0;
};

#endif //MBL_ED_DIAGONALTERM_H
