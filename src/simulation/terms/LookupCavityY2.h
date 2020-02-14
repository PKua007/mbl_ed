//
// Created by Piotr Kubala on 14/02/2020.
//

#ifndef MBL_ED_LOOKUPCAVITYY2_H
#define MBL_ED_LOOKUPCAVITYY2_H

#include <utility>

#include "simulation/CavityConstants.h"
#include "simulation/DoubleHoppingTerm.h"
#include "utils/Assertions.h"

class LookupCavityY2 : public DoubleHoppingTerm {
private:
    double U1{};
    CavityConstants cavityConstants;
    CavityConstants::Realisation currentRealisation;

public:
    LookupCavityY2(double U1, CavityConstants cavityConstants, std::size_t realisationIndex = 0)
            : U1{U1}, cavityConstants{std::move(cavityConstants)}
    {
        Expects(U1 >= 0);
        Expects(realisationIndex < this->cavityConstants.size());
        this->currentRealisation = this->cavityConstants[realisationIndex];
    }

    double calculate(const HopData &firstHopData, const HopData &secondHopData,
                     const HamiltonianGenerator &generator) override;

    void changeRealisation(std::size_t index);
};


#endif //MBL_ED_LOOKUPCAVITYY2_H
