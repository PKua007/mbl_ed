//
// Created by pkua on 02.05.2020.
//

#ifndef MBL_ED_EDEVOLVER_H
#define MBL_ED_EDEVOLVER_H

#include "Evolver.h"

#include "core/Eigensystem.h"

/**
 * @brief Evolver using exact diagonalization technique.
 */
class EDEvolver : public Evolver {
private:
    const Eigensystem &eigensystem;
    arma::cx_mat evolutionOperator;
    arma::cx_vec currentState;
    double dt{};
    double t{};
    std::size_t currentStep{};
    std::size_t numSteps{};

public:
    /**
     * @brief Constructs the evolver using the given @a eigensystem.
     */
    explicit EDEvolver(const Eigensystem &eigensystem);

    void prepareFor(const arma::cx_vec &initialState, double maxTime, std::size_t numSteps_) override;
    void evolve() override;
    [[nodiscard]] const arma::cx_vec &getCurrentState() const override;
    [[nodiscard]] double getDt() const override;
};


#endif //MBL_ED_EDEVOLVER_H
