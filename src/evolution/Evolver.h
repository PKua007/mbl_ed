//
// Created by pkua on 02.05.2020.
//

#ifndef MBL_ED_EVOLVER_H
#define MBL_ED_EVOLVER_H

#include <armadillo>

class Evolver {
public:
    virtual ~Evolver() = default;

    virtual void prepareFor(const arma::cx_vec &initialState, double tMax,
                            std::size_t steps) = 0;
    virtual void evolve() = 0;
    [[nodiscard]] virtual const arma::cx_vec &getCurrentState() const = 0;
    [[nodiscard]] virtual double getCurrentTime() const = 0;
};

#endif //MBL_ED_EVOLVER_H
