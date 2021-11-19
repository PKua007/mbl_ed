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
    enum class StoreLevel {
        NONE,
        EIGENENERGIES,
        EIGENSYSTEM
    };

    /**
     * @brief If true, the diagonalization will produce both eigenvalues and eigenvectors.
     */
    bool calculateEigenvectors{};

    /**
     * @brief How much should be stored to file(s).
     */
    StoreLevel storeLevel = StoreLevel::NONE;

    /**
     * @brief A prefix for all output files produces by the simulation.
     */
    std::string fileSignature{};

};

#endif //MBL_ED_EXACTDIAGONALIZATIONPARAMETERS_H
