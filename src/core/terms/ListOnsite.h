//
// Created by Piotr Kubala on 20/03/2020.
//

#ifndef MBL_ED_LISTONSITE_H
#define MBL_ED_LISTONSITE_H

#include <vector>

#include "core/DiagonalTerm.h"

/**
 * @brief The class representing constant potential in each site from a list provided upfront.
 */
class ListOnsite : public DiagonalTerm {
private:
    std::vector<double> onsitePotential;

public:
    explicit ListOnsite(std::vector<double> onsitePotential) : onsitePotential{std::move(onsitePotential)} { }

    double calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) override;
};



#endif //MBL_ED_LISTONSITE_H
