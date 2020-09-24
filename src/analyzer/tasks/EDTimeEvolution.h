//
// Created by Piotr Kubala on 17/02/2020.
//

#ifndef MBL_ED_EDTIMEEVOLUTION_H
#define MBL_ED_EDTIMEEVOLUTION_H

#include "analyzer/BulkAnalyzerTask.h"
#include "utils/Assertions.h"
#include "evolution/TimeEvolution.h"
#include "evolution/TimeEvolutionParameters.h"

/**
 * @brief BulkAnalyzerTask, which performs time evolution of observables specified by TimeEvolutionParameters from the
 * constructor for a few given Fock states using evolution operator.
 */
class EDTimeEvolution : public BulkAnalyzerTask {
private:
    TimeEvolution timeEvolution;

public:
    /**
     * @brief Creates the class, which will perform evolution for times specified in @a parameters.
     * @details Consult TimeEntryParameters for the description of all parameters. The evolution will be performed
     * separately for all TimeEntryParameters::vectorsToEvolve.
     */
    EDTimeEvolution(const TimeEvolutionParameters &parameters,
                    std::unique_ptr<OservablesTimeEvolution> observablesEvolution);

    /**
     * @brief Adds another Eigensystem to the analyzis, meaning adding another values to means of observables.
     */
    void analyze(const Eigensystem &eigensystem, Logger &logger) override;

    [[nodiscard]] std::string getName() const override;

    /**
     * @brief Stores the result to @a out in the form of a table with a header.
     * @details See TimeEvolution::storeResult for details.
     */
    void storeResult(std::ostream &out) const override;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_EDTIMEEVOLUTION_H
