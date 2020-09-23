//
// Created by Piotr Kubala on 19/02/2020.
//

#ifndef MBL_ED_OCCUPATIONEVOLUTION_H
#define MBL_ED_OCCUPATIONEVOLUTION_H

#include <memory>
#include <utility>

#include "SymmetricMatrix.h"
#include "core/FockBase.h"
#include "Evolver.h"
#include "EvolutionTimeSegment.h"
#include "utils/Logger.h"
#include "PrimaryObservable.h"
#include "SecondaryObservable.h"
#include "CorrelationsTimeEntry.h"

/**
 * @brief The class performing time evolution of @f$ \left< \hat{n}_i \right> @f$ and
 * @f$ \left< \hat{n}_i \hat{n}_j \right> @f$ expected values of observables, where @f$ \hat{n}_i @f$ is onsite number
 * of particles observable.
 */
class OccupationEvolution {
private:
    std::size_t timeStep{};
    double time{};

    std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables;
    std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables;
    std::vector<std::shared_ptr<Observable>> storedObservables;

    std::size_t numOfStoredValues{};

    [[nodiscard]] std::vector<CorrelationsTimeEntry> performTimeSegmentEvolution(std::size_t numSteps, Evolver &evolver,
                                                                                 Logger &logger);

public:
    OccupationEvolution() = default;
    explicit OccupationEvolution(std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables,
                                 std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables,
                                 std::vector<std::shared_ptr<Observable>> storedObservables);
    virtual ~OccupationEvolution() = default;

    /**
     * @brief Perform the evolution of the fock state of index @a initialFockStateIdx for times specified by
     * @a timeSegmentation.
     * @details The actual evolution is performed by a given Evolver. Time segmentation is described in
     * CorrelationsTimeEvolutionParameters.
     * @return The vector of Occupations, where elements corresponds to expected values of observables in subsequent
     * time steps.
     */
    [[nodiscard]] virtual std::vector<CorrelationsTimeEntry>
    perform(const std::vector<EvolutionTimeSegment> &timeSegmentation, const arma::cx_vec &initialState,
            Evolver &evolver, Logger &logger);
};


#endif //MBL_ED_OCCUPATIONEVOLUTION_H
