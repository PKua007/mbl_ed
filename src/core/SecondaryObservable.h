//
// Created by Piotr Kubala on 20/09/2020.
//

#ifndef MBL_ED_SECONDARYOBSERVABLE_H
#define MBL_ED_SECONDARYOBSERVABLE_H

#include <vector>
#include <memory>
#include <type_traits>

#include "Observable.h"
#include "PrimaryObservable.h"

/**
 * @brief An observable, which is calculated from some PrimaryObservable -s.
 */
class SecondaryObservable : public Observable {
public:
    /**
     * @brief Helper method which searches over @a primaryObservables and if found, returns the reference to
     * @a ConcretePrimaryObservable.
     * @details It may be used by derived classes for selecting the observables they need in calculateForObservables().
     */
    template<typename ConcretePrimaryObservable>
    static const ConcretePrimaryObservable &
    findObservable(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) {
        static_assert(std::is_base_of_v<PrimaryObservable, ConcretePrimaryObservable>);
        for (const auto &observable : primaryObservables) {
            try {
                return dynamic_cast<const ConcretePrimaryObservable &>(*observable);
            } catch (std::bad_cast &) { }
        }
        throw std::runtime_error("Desired PrimaryObservable wasn't found");
    }

    virtual void calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) = 0;
};


#endif //MBL_ED_SECONDARYOBSERVABLE_H
