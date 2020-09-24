//
// Created by Piotr Kubala on 20/09/2020.
//

#ifndef MBL_ED_PRIMARYOBSERVABLE_H
#define MBL_ED_PRIMARYOBSERVABLE_H

#include <armadillo>

#include "Observable.h"

/**
 * @brief An observable, which is calculated directly from a state.
 */
class PrimaryObservable : public Observable {
public:
    virtual void calculateForState(const arma::cx_vec &state) = 0;
};


#endif //MBL_ED_PRIMARYOBSERVABLE_H
