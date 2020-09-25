//
// Created by Piotr Kubala on 23/09/2020.
//

#ifndef MBL_ED_OBSERVABLESBUILDER_H
#define MBL_ED_OBSERVABLESBUILDER_H

#include <vector>
#include <string>
#include <memory>

#include "Parameters.h"
#include "core/FockBasis.h"
#include "evolution/PrimaryObservable.h"
#include "evolution/SecondaryObservable.h"

/**
 * @brief A class building PrimaryObservables, SecondaryObservables and stored Observable -s lists from string
 * representations.
 * @details See the code to see what is possible
 */
class ObservablesBuilder {
private:
    std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables;
    std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables;
    std::vector<std::shared_ptr<Observable>> storedObservables;

public:
    /**
     * @brief Bulids the lists of observables.
     * @details All observables represented by @a observables are of course stored Observables, irregardless if they
     * are primary or secondary. Only the minimal set of required PrimaryObservables is creating if not specified
     * explicitly, but needed by SecondaryObservables.
     */
    void build(const std::vector<std::string> &observables, const Parameters &params,
               const std::shared_ptr<FockBasis> &fockBasis);

    /**
     * @brief Returns the list of PrimaryObservables and releases it from the memory of the class.
     */
    std::vector<std::shared_ptr<PrimaryObservable>> releasePrimaryObservables() {
        return std::move(this->primaryObservables);
    }

    /**
     * @brief Returns the list of SecondaryObservable and releases it from the memory of the class.
     */
    std::vector<std::shared_ptr<SecondaryObservable>> releaseSecondaryObservables() {
        return std::move(this->secondaryObservables);
    }

    /**
     * @brief Returns the list of stored Observables and releases it from the memory of the class.
     */
    std::vector<std::shared_ptr<Observable>> releaseStoredObservables() {
        return std::move(this->storedObservables);
    }
};


#endif //MBL_ED_OBSERVABLESBUILDER_H
