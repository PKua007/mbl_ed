//
// Created by Piotr Kubala on 27/01/2020.
//

#ifndef MBL_ED_EXACTDIAGONALIZATIONPARAMETERS_H
#define MBL_ED_EXACTDIAGONALIZATIONPARAMETERS_H

#include "SimulationsSpan.h"

/**
 * @brief The paramaters of Simulation.
 */
struct ExactDiagonalizationParameters {
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

#endif //MBL_ED_EXACTDIAGONALIZATIONPARAMETERS_H
