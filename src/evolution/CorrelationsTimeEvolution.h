//
// Created by pkua on 08.05.2020.
//

#ifndef MBL_ED_CORRELATIONSTIMEEVOLUTION_H
#define MBL_ED_CORRELATIONSTIMEEVOLUTION_H

#include <memory>

#include "SymmetricMatrix.h"
#include "OccupationEvolution.h"
#include "CorrelationsTimeEntry.h"
#include "Evolver.h"
#include "CorrelationsTimeEvolutionParameters.h"

class CorrelationsTimeEvolution {
private:
    /**
     * @brief Evolution of observables for the specifit initial Fock state.
     */
    struct VectorEvolution {
        FockBase::Vector initialVector;
        std::vector<CorrelationsTimeEntry> timeEntries{};

        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string getInitialVectorSignature() const;
    };

    std::shared_ptr<FockBase> fockBase;
    std::size_t marginSize{};
    std::size_t numberOfSites{};
    std::vector<VectorEvolution> vectorEvolutions{};
    double maxTime{};
    std::size_t numSteps{};

    [[nodiscard]] std::size_t getNumberOfSites() const;

public:
    /**
     * @brief Creates the class, which will perform evolution from 0 to @a maxTime time, dividing it into @a numSteps
     * steps.
     * @details Consult CorrelationTimeEntry for the description of @a marginSize. The evolution will be performed
     * separately for all @a vectorsToEvolve. To know what observables are measures, check CorrelationTimeEntry.
     */
    explicit CorrelationsTimeEvolution(const CorrelationsTimeEvolutionParameters &parameters);

    /**
     * @brief Adds another Eigensystem to the analyzis, to all observables means are enriched by more data.
     * @details The number of sites and number of particles are determined in the first invocation and must be kept the
     * same in next ones.
     */
    void addEvolution(Evolver &evolver, std::ostream &logger);


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
    void storeResult(std::ostream &out) const;
};


#endif //MBL_ED_CORRELATIONSTIMEEVOLUTION_H
