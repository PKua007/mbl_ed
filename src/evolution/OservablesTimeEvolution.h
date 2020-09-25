//
// Created by Piotr Kubala on 19/02/2020.
//

#ifndef MBL_ED_OSERVABLESTIMEEVOLUTION_H
#define MBL_ED_OSERVABLESTIMEEVOLUTION_H

#include <memory>
#include <utility>

#include "SymmetricMatrix.h"
#include "core/FockBasis.h"
#include "Evolver.h"
#include "EvolutionTimeSegment.h"
#include "utils/Logger.h"
#include "PrimaryObservable.h"
#include "SecondaryObservable.h"
#include "TimeEvolutionEntry.h"

/**
 * @brief The class performing time evolution of all PrimaryObservable -s and SecondaryObservable -s specified by the
 * setters.
 */
class OservablesTimeEvolution {
private:
    std::size_t timeStep{};
    double time{};

    std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables;
    std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables;
    std::vector<std::shared_ptr<Observable>> storedObservables;

    std::size_t numOfObservableValues{};

    [[nodiscard]] std::vector<TimeEvolutionEntry> performTimeSegmentEvolution(std::size_t numSteps, Evolver &evolver,
                                                                              Logger &logger);

public:
    virtual ~OservablesTimeEvolution() = default;

    /**
     * @brief Sets which PrimaryObservable -s will be calculated.
     */
    void setPrimaryObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables_) {
        this->primaryObservables = primaryObservables_;
    }

    /**
     * @brief Sets which SecondaryObservable -s will be calculated.
     */
    void setSecondaryObservables(const std::vector<std::shared_ptr<SecondaryObservable>> &secondaryObservables_) {
        this->secondaryObservables = secondaryObservables_;
    }

    /**
     * @brief Sets which Observable -s will be return as a result of perform() invokation.
     */
    void setStoredObservables(const std::vector<std::shared_ptr<Observable>> &storedObservables_);

    /**
     * @brief Perform the evolution of the @a initialFockStateIdx for times specified by @a timeSegmentation.
     * @details The actual evolution is performed by the given Evolver. Time segmentation is described in
     * TimeEvolutionParameters.
     * @return The vector of TimeEvolutionEntry -ies, where elements corresponds to expected values of observables in
     * subsequent time steps. The observables to be returned are determined by setStoredObservables().
     */
    [[nodiscard]] virtual std::vector<TimeEvolutionEntry>
    perform(const std::vector<EvolutionTimeSegment> &timeSegmentation, const arma::cx_vec &initialState,
            Evolver &evolver, Logger &logger);
};


#endif //MBL_ED_OSERVABLESTIMEEVOLUTION_H
