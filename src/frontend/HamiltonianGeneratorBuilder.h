//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_HAMILTONIANGENERATORBUILDER_H
#define MBL_ED_HAMILTONIANGENERATORBUILDER_H

#include <memory>

#include "simulation/HamiltonianGenerator.h"
#include "simulation/RND.h"
#include "Parameters.h"


class HamiltonianGeneratorBuilder {
public:
    std::unique_ptr<HamiltonianGenerator> build(const Parameters &params, std::shared_ptr<FockBase> fockBase, RND &rnd);
};


#endif //MBL_ED_HAMILTONIANGENERATORBUILDER_H
