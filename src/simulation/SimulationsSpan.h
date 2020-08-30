//
// Created by Piotr Kubala on 20/08/2020.
//

#ifndef MBL_ED_SIMULATIONSSPAN_H
#define MBL_ED_SIMULATIONSSPAN_H

#include <cstddef>

struct SimulationsSpan {
    /**
     * @brief The starting index (inclusive) of the simulation (from 0).
     */
    std::size_t from{};

    /**
     * @brief The final index (exclusive) of the simulation.
     */
    std::size_t to{};

    /**
     * @brief The total number of simulations. SimulationParameters::to and SimulationParameters::from parameters
     * describe the actual range of simulations to be performed in this run, while total number refers to all runs.
     */
    std::size_t total{};
};


#endif //MBL_ED_SIMULATIONSSPAN_H
