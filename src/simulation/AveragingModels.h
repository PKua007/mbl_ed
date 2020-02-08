//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_AVERAGINGMODELS_H
#define MBL_ED_AVERAGINGMODELS_H

#include <cmath>

#include "utils/Assertions.h"
#include "RND.h"

/**
 * @brief Averaging model, used in Simulation, which samples new onsite disorder terms for each simulation.
 * @tparam HamiltonianGenerator_t concrete hamiltonian generator that it to be prepared for new average entry
 */
template<typename HamiltonianGenerator_t>
class OnsiteDisorderAveragingModel {
public:
    static void setupHamiltonianGenerator(HamiltonianGenerator_t &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        static_cast<void>(simulationIndex);
        static_cast<void>(numberOfSimulations);
        hamiltonianGenerator.resampleOnsiteEnergies(rnd);
    }
};

/**
 * @brief Averaging model, used in Simulation, which distributes evenly phi0 values across simulations based on
 * @a simulationIndex and @a numberOfSimulations.
 * @tparam HamiltonianGenerator_t concrete hamiltonian generator that it to be prepared for new average entry
 */
template<typename HamiltonianGenerator_t>
class UniformPhi0AveragingModel {
public:
    static void setupHamiltonianGenerator(HamiltonianGenerator_t &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        Expects(numberOfSimulations > 0);
        Expects(simulationIndex < numberOfSimulations);

        double phi0 = 2*M_PI*simulationIndex/numberOfSimulations;
        hamiltonianGenerator.resampleOnsiteEnergies(rnd);
        hamiltonianGenerator.setPhi0(phi0);
    }
};

/**
 * @brief Averaging model, used in Simulation, which picks random phi0 values across simulations based on
 * @a simulationIndex and @a numberOfSimulations.
 * @tparam HamiltonianGenerator_t concrete hamiltonian generator that it to be prepared for new average entry
 */
template<typename HamiltonianGenerator_t>
class RandomPhi0AveragingModel {
public:
    static void setupHamiltonianGenerator(HamiltonianGenerator_t &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        Expects(numberOfSimulations > 0);
        Expects(simulationIndex < numberOfSimulations);

        hamiltonianGenerator.resampleOnsiteEnergies();
        hamiltonianGenerator.setPhi0(2*M_PI*rnd());
    }
};

#endif //MBL_ED_AVERAGINGMODELS_H
