//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_HOPPINGTERM_H
#define MBL_ED_HOPPINGTERM_H

#include "FockBase.h"

class HamiltonianGenerator;

class HoppingTerm {
public:
    virtual ~HoppingTerm() = default;

    virtual double calculate(const FockBase::Vector &from, const FockBase::Vector &to, std::size_t fromSite,
                             std::size_t toSite, const HamiltonianGenerator &generator) = 0;
};

#endif //MBL_ED_HOPPINGTERM_H
