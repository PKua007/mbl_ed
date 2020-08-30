//
// Created by Piotr Kubala on 09/02/2020.
//

#include "HubbardHop.h"
#include "core/HamiltonianGenerator.h"

double HubbardHop::calculate(const HopData &hopData, const HamiltonianGenerator &generator)
{
    Expects(generator.getSiteDistance(hopData.fromSite, hopData.toSite) == 1);
    static_cast<void>(hopData.fromVector);
    static_cast<void>(hopData.toVector);

    return -this->J;
}

HubbardHop::HubbardHop(double J) : J{J} {
    Expects(J >= 0);
}
