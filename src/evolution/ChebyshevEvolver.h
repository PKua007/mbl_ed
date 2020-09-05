//
// Created by pkua on 02.05.2020.
//

#ifndef MBL_ED_CHEBYSHEVEVOLVER_H
#define MBL_ED_CHEBYSHEVEVOLVER_H


#include "Evolver.h"

/**
 * @brief Evolver usign Chebyshev expansion technique from paper:
 * <em>Many-body localization in presence of cavity mediated long-range interactions</em>
 */
class ChebyshevEvolver : public Evolver {
private:
    const arma::sp_mat &hamiltonian;
    arma::cx_vec currentState;
    double a{};
    double b{};
    std::size_t N{};
    double t{};
    double dt{};
    std::size_t currentStep{};
    std::size_t maxSteps{};
    std::ostream &logger;

    static constexpr std::size_t MIN_EIGVAL = 20;
    static constexpr double MAXIMAL_NORM_LEAKAGE = 1e-12;

    void findSpectrumRange();
    void optimizeOrder(const arma::cx_vec &initialState);
    [[nodiscard]] arma::cx_vec evolveState(const arma::cx_vec &state);

public:
    /**
     * @brief Constructs the evolver which will be using given @a hamiltonian
     */
    ChebyshevEvolver(const arma::sp_mat &hamiltonian, std::ostream &logger);

    void prepareFor(const arma::cx_vec &initialState, double maxTime, std::size_t maxSteps_) override;
    void evolve() override;
    [[nodiscard]] const arma::cx_vec &getCurrentState() const override;
    [[nodiscard]] double getDt() const override;
};


#endif //MBL_ED_CHEBYSHEVEVOLVER_H
