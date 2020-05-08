//
// Created by Piotr Kubala on 17/02/2020.
//

#ifndef MBL_ED_EDCORRELATIONSTIMEEVOLUTION_H
#define MBL_ED_EDCORRELATIONSTIMEEVOLUTION_H

#include "analyzer/BulkAnalyzerTask.h"
#include "utils/Assertions.h"
#include "evolution/CorrelationsTimeEvolution.h"

/**
 * @brief BulkAnalyzerTask, which performs time evolution of correlations and fluctuations for a few given Fock states.
 * @details Consult CorrelationsTimeEntry and its helper classes for the description.
 */
class EDCorrelationsTimeEvolution : public BulkAnalyzerTask {
private:
    CorrelationsTimeEvolution correlationsTimeEvolution;

public:
    /**
     * @brief Creates the class, which will perform evolution from 0 to @a maxTime time, dividing it into @a numSteps
     * steps.
     * @details Consult CorrelationTimeEntry for the description of @a marginSize. The evolution will be performed
     * separately for all @a vectorsToEvolve. To know what observables are measures, check CorrelationTimeEntry.
     */
    EDCorrelationsTimeEvolution(double maxTime, std::size_t numSteps, std::size_t numbefOfSites, std::size_t marginSize,
                                std::shared_ptr<FockBase> fockBase, const std::vector<FockBase::Vector> &vectorsToEvolve);

    /**
     * @brief Adds another Eigensystem to the analyzis, to all observables means are enriched by more data.
     * @details The number of sites and number of particles are determined in the first invocation and must be kept the
     * same in next ones.
     */
    void analyze(const Eigensystem &eigensystem, std::ostream &logger) override;

    [[nodiscard]] std::string getName() const override;

    /**
     * @brief Stores the result to @a out in the form of a table with a header.
     * @details It is constructed from horizontally glued tables for all vectors to evolve specified in the constructor.
     * <p> The table for a single vector is constructed from CorrelationsTimeEntry -ies. The header is obtained from
     * CorrelationsTimeEntry::getHeader and the values are for subsequent times and are constructed by
     * CorrelationsTimeEntry::toString. Moreover, to the first column name, in header, namely time, a representation
     * of the initial vector is prepended in the form of site1occupaiton.site2occupaiton.[...].lastSiteOccupation_.
     * <p> To sum up, for 2 vectors, {1, 1, 1, 1} and {2, 0, 2, 0}, for @a maxTime = 1 and @a numSteps = 3, it will
     * look like:
     * <pre>
     * 1.1.1.1_t [observables] 2.0.2.0_t [observables]
     * 0 [values] 0 [values]
     * 0.5 [values] 0.5 [values]
     * 1 [values] 1 [values]
     * </pre>
     */
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_EDCORRELATIONSTIMEEVOLUTION_H
