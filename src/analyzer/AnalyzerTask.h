//
// Created by Piotr Kubala on 16/01/2020.
//

#ifndef MBL_ED_ANALYZERTASK_H
#define MBL_ED_ANALYZERTASK_H

#include <vector>
#include <string>

#include "core/Eigensystem.h"
#include "simulation/Restorable.h"
#include "utils/Logger.h"

/**
 * @brief An analyzis task, which will be performed on each Eigensystem given to Analyzer.
 */
class AnalyzerTask : public Restorable {
public:
    virtual ~AnalyzerTask() = default;

    /**
     * @brief Performs the analyzis on a given @a eigensystem.
     *
     * Note, that it should be able to accept multiple Eigensystem -s from subsequent simulations. Whether it will
     * average or just append new results it is up to a specific AnalyzerTask.
     */
    virtual void analyze(const Eigensystem &eigensystem, Logger &logger) = 0;

    /**
     * @brief Returns the name of the analyzer task. It can be used for example for file name suffixes.
     */
    [[nodiscard]] virtual std::string getName() const = 0;
};


#endif //MBL_ED_ANALYZERTASK_H
