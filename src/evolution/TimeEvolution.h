//
// Created by pkua on 08.05.2020.
//

#ifndef MBL_ED_TIMEEVOLUTION_H
#define MBL_ED_TIMEEVOLUTION_H

#include <memory>
#include <variant>

#include "SymmetricMatrix.h"
#include "OservablesTimeEvolution.h"
#include "Evolver.h"
#include "TimeEvolutionParameters.h"
#include "simulation/Restorable.h"
#include "utils/Logger.h"

/**
 * @brief The class supervising the whole process of performing the evolution and calculating the observables for
 * some discrete time steps.
 * @details The subsequents evolutions are performed when addEvolution() is invoked and the overall results are
 * averaged. The evolution is performed for one or more initial states. The observables measured are determined by
 * TimeEvolutionParameters from the constructor.
 */
class TimeEvolution : public Restorable {
private:
    /**
     * @brief Evolution of observables for the specific initial Fock state.
     */
    struct VectorEvolution : public Restorable {
        using ExternalVector = TimeEvolutionParameters::ExternalVector;

        std::variant<FockBasis::Vector, ExternalVector> initialVector;
        std::vector<TimeEvolutionEntry> timeEntries{};
        std::string observablesHeader;

        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string getInitialVectorName() const;
        void storeState(std::ostream &binaryOut) const override;
        void joinRestoredState(std::istream &binaryIn) override;
        void clear() override;
    };

    std::shared_ptr<FockBasis> fockBasis;
    std::unique_ptr<OservablesTimeEvolution> occupationEvolution;
    std::vector<VectorEvolution> vectorEvolutions{};
    std::vector<EvolutionTimeSegment> timeSegmentation{};

public:
    /**
     * @brief Creates the class, which will perform evolution parametrised by @a parameters
     * @details The evolution will be performed separately for all TimeEvolutionParameters::vectorsToEvolve.
     * @a occupationEvolution is used for the acutal evolution process.
     */
    explicit TimeEvolution(const TimeEvolutionParameters &parameters,
                           std::unique_ptr<OservablesTimeEvolution> occupationEvolution);

    /**
     * @brief Adds another evolution to the analyzis.
     * @details The actual evolution is done by OservablesTimeEvolution using the given @a evolver. New data are
     * averaged will the old ones. If field @a initialVectors from @a parameters from the constructor contained some
     * TimeEvolutionParameters::ExternalVector alternatives, the actual arma::cx_vec vectors should be passed through
     * @a externalVectors. The rest are FockBasis::Vectors product vectors and are prepared on the go.
     */
    void addEvolution(Evolver &evolver, Logger &logger, const std::vector<arma::cx_vec> &externalVectors = {});

    /**
     * @brief Stores the result to @a out in the form of a table with a header.
     * @details It is constructed from horizontally glued tables for all vectors to evolve specified in the constructor.
     * <p> The table for a single vector is constructed from TimeEvolutionEntry -ies. The header contains time with
     * evolved vector, in format (site1occupaiton).(site2occupaiton).[...].(lastSiteOccupation)_t, see below, and all
     * stored observables names.
     * <p> Example: for 2 vectors, {1, 1, 1, 1} and {2, 0, 2, 0}, for @a maxTime = 1 and @a numSteps = 3 (see
     * TimeEvolutionParameters), it will look like:
     * <pre>
     * 1.1.1.1_t [observables] 2.0.2.0_t [observables]
     * 0 [values] 0 [values]
     * 0.5 [values] 0.5 [values]
     * 1 [values] 1 [values]
     * </pre>
     */
    void storeResult(std::ostream &out) const;

    /**
     * @brief Returns how many external vector are present in this evolution.
     */
    [[nodiscard]] std::size_t countExternalVectors() const;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_TIMEEVOLUTION_H
