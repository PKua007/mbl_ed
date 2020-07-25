//
// Created by pkua on 09.05.2020.
//

#ifndef MBL_ED_CORRELATIONSTIMEEVOLUTIONPARAMETERS_H
#define MBL_ED_CORRELATIONSTIMEEVOLUTIONPARAMETERS_H

#include <memory>
#include <vector>

#include "simulation/FockBase.h"
#include "utils/Assertions.h"
#include "EvolutionTimeSegment.h"

/**
 * @brief Parameters of time evolution.
 */
struct CorrelationsTimeEvolutionParameters {
    /**
     * @brief How the times are sampled.
     * @details Fox example, if it is {{maxTime = 1, numSteps = 2}, {maxTime = 3, numSteps = 1}}, the resulting times
     * should be 0, 0.5, 1, 3
     */
    std::vector<EvolutionTimeSegment> timeSegmentation;
    std::size_t numberOfSites{};

    /**
     * @brief Margin size for borderless observables (see CorrelationsTimeEntry)
     */
    std::size_t marginSize{};
    std::shared_ptr<FockBase> fockBase{};
    std::vector<FockBase::Vector> vectorsToEvolve{};

    /**
     * @brief Constructs CorrelationsTimeEvolutionParameters::vectorsToEvolve from a string.
     * @details It can be either a tag
     * <ul>
     * <li> unif - 1.1.1.1
     * <li> dw - 2.0.2.0
     * </ul>
     * or occupation representation 2.3.0.0.1
     */
    void setVectorsToEvolveFromTags(const std::vector<std::string> &strings) {
        Expects(this->numberOfSites > 0);

        this->vectorsToEvolve.clear();
        for (const auto &string : strings) {
            try {
                // Try a tag representation
                this->vectorsToEvolve.emplace_back(this->numberOfSites, string);
            } catch (FockVectorParseException &) {
                // If tag failed, try occupation representation
                this->vectorsToEvolve.emplace_back(string);
            }
        }
    }
};

#endif //MBL_ED_CORRELATIONSTIMEEVOLUTIONPARAMETERS_H
