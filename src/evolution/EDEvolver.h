//
// Created by pkua on 02.05.2020.
//

#ifndef MBL_ED_EDEVOLVER_H
#define MBL_ED_EDEVOLVER_H


#include "Evolver.h"

class EDEvolver : public Evolver {
private:
    arma::cx_mat evolutionOperator;
    arma::cx_vec currentState;
    double dt{};
    double t{};

public:
    void prepareFor(const arma::sp_mat &hamiltonian, const arma::cx_vec &initialState, double tMax,
                    std::size_t steps) override;
    void evolve() override;
    [[nodiscard]] const arma::cx_vec &getCurrentState() const override;
    [[nodiscard]] double getCurrentTime() const override;
};


#endif //MBL_ED_EDEVOLVER_H
