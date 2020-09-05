//
// Created by Piotr Kubala on 05/09/2020.
//

#include "HamiltonianGeneratorMother.h"

#include "core/terms/HubbardOnsite.h"
#include "core/terms/HubbardHop.h"
#include "core/terms/QuasiperiodicDisorder.h"

std::unique_ptr<HamiltonianGenerator>
HamiltonianGeneratorMother::hubbardQuasiperiodic(double J, double U, double W, double beta, double phi0, bool usePBC) {
    auto hamiltonianGenerator = std::make_unique<HamiltonianGenerator>(fockBase, usePBC);
    hamiltonianGenerator->addHoppingTerm(std::make_unique<HubbardHop>(J));
    hamiltonianGenerator->addDiagonalTerm(std::make_unique<HubbardOnsite>(U));
    hamiltonianGenerator->addDiagonalTerm(std::make_unique<QuasiperiodicDisorder>(W, beta, phi0));
    return hamiltonianGenerator;
}
