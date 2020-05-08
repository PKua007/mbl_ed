//
// Created by pkua on 02.05.2020.
//

#ifndef MBL_ED_CHEBYSHEVEVOLVER_H
#define MBL_ED_CHEBYSHEVEVOLVER_H


#include "Evolver.h"

class ChebyshevEvolver : public Evolver {
private:
    const arma::sp_mat &hamiltonian;
    std::vector<arma::cx_vec> chebyshevVectors;
    arma::cx_vec currentState;
    double a{};
    double b{};
    std::size_t N{};
    double t{};
    double dt{};
    double Nfactor{};

public:
    explicit ChebyshevEvolver(const arma::sp_mat &hamiltonian, double Nfactor = 1.5);

    void prepareFor(const arma::cx_vec &initialState, double tMax,
                    std::size_t steps) override;
    void evolve() override;
    [[nodiscard]] const arma::cx_vec &getCurrentState() const override;
    [[nodiscard]] double getCurrentTime() const override;
};


#endif //MBL_ED_CHEBYSHEVEVOLVER_H
