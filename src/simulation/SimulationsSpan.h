//
// Created by Piotr Kubala on 20/08/2020.
//

#ifndef MBL_ED_SIMULATIONSSPAN_H
#define MBL_ED_SIMULATIONSSPAN_H

#include <cstddef>

struct SimulationsSpan {
    std::size_t from{};
    std::size_t to{};
    std::size_t totalSimulations{};
};


#endif //MBL_ED_SIMULATIONSSPAN_H
