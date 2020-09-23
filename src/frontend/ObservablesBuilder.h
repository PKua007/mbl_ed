//
// Created by Piotr Kubala on 23/09/2020.
//

#ifndef MBL_ED_OBSERVABLESBUILDER_H
#define MBL_ED_OBSERVABLESBUILDER_H

#include <vector>
#include <string>
#include <memory>

#include "Parameters.h"
#include "core/FockBase.h"
#include "evolution/PrimaryObservable.h"
#include "evolution/SecondaryObservable.h"


class ObservablesBuilder {
private:
    std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables;
    std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables;
    std::vector<std::shared_ptr<Observable>> storedObservables;

public:
    void build(const std::vector<std::string> &observables, const Parameters &params,
               const std::shared_ptr<FockBase> &fockBase);

    std::vector<std::shared_ptr<PrimaryObservable>> releasePrimaryObservables() {
        return std::move(this->primaryObservables);
    }

    std::vector<std::shared_ptr<SecondaryObservable>> releaseSecondaryObservables() {
        return std::move(this->secondaryObservables);
    }

    std::vector<std::shared_ptr<Observable>> releaseStoredObservables() {
        return std::move(this->storedObservables);
    }
};


#endif //MBL_ED_OBSERVABLESBUILDER_H
