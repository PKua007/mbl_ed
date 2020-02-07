//
// Created by Piotr Kubala on 27/01/2020.
//

#ifndef MBL_ED_SIMULATIONPARAMETERS_H
#define MBL_ED_SIMULATIONPARAMETERS_H

/**
 * @brief The paramaters of Simulation.
 */
struct SimulationParameters {
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
    std::size_t totalSimulations{};

    /**
     * @brief If true, the diagonalization will produce both eigenvalues and eigenvectors.
     */
    bool calculateEigenvectors{};

    /**
     * @brief If true, after each diagonalization eigenenergies will be saved to a file.
     */
    bool saveEigenenergies{};

    /**
     * @brief A prefix for all output files produces by the simulation.
     */
    std::string fileSignature{};
};

#endif //MBL_ED_SIMULATIONPARAMETERS_H
