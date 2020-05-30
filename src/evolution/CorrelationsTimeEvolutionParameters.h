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
     * @brief Constructs CorrelationsTimeEvolutionParameters::vectorsToEvolve from a string @a tag.
     * @details Supported tags (for 4 sites), giving appropriate vectors:
     * <ul>
     * <li> unif - 1.1.1.1
     * <li> dw - 2.0.2.0
     * <li> both - 1.1.1.1, 2.0.2.0
     * </ul>
     */
    void setVectorsToEvolveFromTag(const std::string &tag) {
        Expects(this->numberOfSites > 0);

        FockBase::Vector uniform(this->numberOfSites, 1);
        FockBase::Vector densityWave(this->numberOfSites);
        for (std::size_t i{}; i < densityWave.size(); i += 2)
            densityWave[i] = 2;

        if (tag == "unif") {
            this->vectorsToEvolve = {uniform};
        } else if (tag == "dw") {
            this->vectorsToEvolve = {densityWave};
        } else if (tag == "both") {
            this->vectorsToEvolve = {uniform, densityWave};
        } else {
            throw std::runtime_error("vectors to evolve must be unif/dw/both"
                                     "\nunif - 1.1.1.1; dw - 2.0.2.0; both - both ;)");
        }
    }
};

#endif //MBL_ED_CORRELATIONSTIMEEVOLUTIONPARAMETERS_H
