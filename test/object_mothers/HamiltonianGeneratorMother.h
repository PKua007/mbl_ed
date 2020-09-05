//
// Created by Piotr Kubala on 05/09/2020.
//

#ifndef MBL_ED_HAMILTONIANGENERATORMOTHER_H
#define MBL_ED_HAMILTONIANGENERATORMOTHER_H

#include <memory>

#include "core/HamiltonianGenerator.h"
#include "core/FockBaseGenerator.h"


class HamiltonianGeneratorMother {
private:
    std::shared_ptr<FockBase> fockBase;

public:
    HamiltonianGeneratorMother(std::size_t N, std::size_t K) : fockBase{FockBaseGenerator{}.generate(N, K)}
    { }

   std::unique_ptr<HamiltonianGenerator> hubbardQuasiperiodic(double J, double U, double W, double beta, double phi0,
                                              bool usePBC = false);
};


#endif //MBL_ED_HAMILTONIANGENERATORMOTHER_H
