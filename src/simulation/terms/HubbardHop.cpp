//
// Created by Piotr Kubala on 09/02/2020.
//

#include "HubbardHop.h"
#include "simulation/HamiltonianGenerator.h"

double HubbardHop::calculate(const FockBase::Vector &from, const FockBase::Vector &to, std::size_t fromSite,
                             std::size_t toSite, const HamiltonianGenerator &generator)
{
    Expects(generator.getSiteDistance(fromSite, toSite) == 1);
    static_cast<void>(from);
    static_cast<void>(to);

    return -this->J;
}

HubbardHop::HubbardHop(double J) : J{J} {
    Expects(J >= 0);
}
