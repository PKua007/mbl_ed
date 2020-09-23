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

class SecondaryObservable : public Observable {
public:
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
