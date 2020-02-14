//
// Created by Piotr Kubala on 14/02/2020.
//

#ifndef MBL_ED_DOUBLEHOPPINGTERM_H
#define MBL_ED_DOUBLEHOPPINGTERM_H

class HamiltonianGenerator;
class HopData;

class DoubleHoppingTerm {
public:
    virtual ~DoubleHoppingTerm() = default;

    virtual double calculate(const HopData &firstHopData, const HopData &secondHopData,
                             const HamiltonianGenerator &generator) = 0;
};

#endif //MBL_ED_DOUBLEHOPPINGTERM_H
