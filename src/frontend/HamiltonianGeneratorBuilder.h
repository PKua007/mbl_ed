//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_HAMILTONIANGENERATORBUILDER_H
#define MBL_ED_HAMILTONIANGENERATORBUILDER_H

#include <memory>

#include "core/HamiltonianGenerator.h"
#include "core/RND.h"
#include "Parameters.h"


/**
 * @brief A class responsible for building HamiltonianGenerator based on Parameters.
 * @details It is extracted from Frontend mainly for decreasing compilation times.
 */
class HamiltonianGeneratorBuilder {
public:
    std::unique_ptr<HamiltonianGenerator> build(const Parameters &params, std::shared_ptr<FockBase> fockBase, RND &rnd);
};


#endif //MBL_ED_HAMILTONIANGENERATORBUILDER_H
