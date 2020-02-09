//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_HUBBARDHOP_H
#define MBL_ED_HUBBARDHOP_H


#include "utils/Assertions.h"
#include "simulation/HoppingTerm.h"

class HubbardHop : public HoppingTerm {
private:
    double J{};

public:
    explicit HubbardHop(double J);

    /**
     * @brief Returns constant -J (see CavityHamiltonianParameters) for a hop between neighbouring
     * sites (according to PBC or OBC) and 0 for larger hops.
     */
    double calculate(const FockBase::Vector &from, const FockBase::Vector &to, std::size_t fromSite,
                     std::size_t toSite, const HamiltonianGenerator &generator) override;
};

#endif //MBL_ED_HUBBARDHOP_H
