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
    double maxTimeForN{};
    double t{};
    double dt{};
    double Nfactor{};
    std::size_t step{};
    std::size_t steps{};

    void rebuildChebychevVectors(const arma::cx_vec &initialState);

public:
    explicit ChebyshevEvolver(const arma::sp_mat &hamiltonian, std::size_t N = 1000, double Nfactor = 1.5);

    void prepareFor(const arma::cx_vec &initialState, double maxTime,
                    std::size_t steps) override;
    void evolve() override;
    [[nodiscard]] const arma::cx_vec &getCurrentState() const override;

    [[nodiscard]] double getMaxTimeForN() const;
};


#endif //MBL_ED_CHEBYSHEVEVOLVER_H
