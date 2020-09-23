//
// Created by Piotr Kubala on 20/09/2020.
//

#ifndef MBL_ED_SECONDARYOBSERVABLE_H
#define MBL_ED_SECONDARYOBSERVABLE_H

#include <vector>
#include <memory>

#include "Observable.h"
#include "PrimaryObservable.h"

class SecondaryObservable : public Observable {
public:
    virtual void calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &rawObservables) = 0;
};


#endif //MBL_ED_SECONDARYOBSERVABLE_H
