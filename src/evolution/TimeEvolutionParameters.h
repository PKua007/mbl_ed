//
// Created by pkua on 09.05.2020.
//

#ifndef MBL_ED_TIMEEVOLUTIONPARAMETERS_H
#define MBL_ED_TIMEEVOLUTIONPARAMETERS_H

#include <memory>
#include <vector>
#include <variant>

#include "core/FockBase.h"
#include "EvolutionTimeSegment.h"
#include "PrimaryObservable.h"
#include "SecondaryObservable.h"

/**
 * @brief Parameters of time evolution.
 */
struct TimeEvolutionParameters {
    /**
     * @brief Helper struct representing "an external initial vector", which is not FockBase::Vector and should be
     * passed right before the evolution (see CorrelationsTimeEvolution::addEvolution)
     */
    struct ExternalVector {
        std::string name;
    };

    /**
     * @brief How the times are sampled.
     * @details Fox example, if it is {{maxTime = 1, numSteps = 2}, {maxTime = 3, numSteps = 1}}, the resulting times
     * should be 0, 0.5, 1, 3
     */
    std::vector<EvolutionTimeSegment> timeSegmentation;

    std::size_t numberOfSites{};
    std::shared_ptr<FockBase> fockBase{};

    /**
     * @brief Initial vectors to evolve.
     * @details Initial vector can be either a Fock basis vector known from the beginning or "an external vector",
     * which is not known from the beginning and should be passed before evolution (see
     * CorrelationsTimeEvolution::addEvolution)
     */
    std::vector<std::variant<FockBase::Vector, ExternalVector>> vectorsToEvolve{};

    std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables;
    std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables;
    std::vector<std::shared_ptr<Observable>> storedObservables;

    /**
     * @brief Constructs CorrelationsTimeEvolutionParameters::vectorsToEvolve from a string.
     * @details It can be either a tag
     * <ul>
     * <li> unif - 1.1.1.1
     * <li> dw - 2.0.2.0
     * </ul>
     * or occupation representation 2.3.0.0.1
     */
    void setVectorsToEvolveFromTags(const std::vector<std::string> &strings);

    [[nodiscard]] std::size_t countStoredObservableValues() const;
    [[nodiscard]] std::string generateStoredObservablesHeader() const;
};

#endif //MBL_ED_TIMEEVOLUTIONPARAMETERS_H
