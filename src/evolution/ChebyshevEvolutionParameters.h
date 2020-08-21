//
// Created by Piotr Kubala on 20/08/2020.
//

#ifndef MBL_ED_CHEBYSHEVEVOLUTIONPARAMETERS_H
#define MBL_ED_CHEBYSHEVEVOLUTIONPARAMETERS_H

#include <string>

#include "CorrelationsTimeEvolutionParameters.h"

struct ChebyshevEvolutionParameters {
    std::size_t from{};
    std::size_t to{};
    std::size_t totalSimulations{};
    CorrelationsTimeEvolutionParameters timeEvolutionParameters;
};


#endif //MBL_ED_CHEBYSHEVEVOLUTIONPARAMETERS_H
