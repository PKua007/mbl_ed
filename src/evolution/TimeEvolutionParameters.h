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
     * passed right before the evolution (see TimeEvolution::addEvolution)
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

    /**
     * @brief PrimaryObservable -s which should be calculated during the evolution.
     */
    std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables;

    /**
     * @brief SecondaryObservable -s which should be calculated during the evolution.
     */
    std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables;

    /**
     * @brief Observables (may be both PrimaryObservables and SecondaryObservables) which should be in the output
     * of the evolution.
     */
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

    /**
     * @brief Returns how many observable values (so summing sets from all stored observables) will be in the output.
     */
    [[nodiscard]] std::size_t countStoredObservableValues() const;

    /**
     * @brief Constructs space-separated header of stored observable names from all stored observables sets,
     * preserving the order.
     * @brief Namely, it will look like:
     * <pre>
     * [observable 1 value 1] ... [o. 1 last value] ... [last observable value 1] ... [l. o. last value]
     * </pre>
     */
    [[nodiscard]] std::string generateStoredObservablesHeader() const;
};

#endif //MBL_ED_TIMEEVOLUTIONPARAMETERS_H
