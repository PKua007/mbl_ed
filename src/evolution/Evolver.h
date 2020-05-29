//
// Created by pkua on 02.05.2020.
//

#ifndef MBL_ED_EVOLVER_H
#define MBL_ED_EVOLVER_H

#include <armadillo>

/**
 * @brief A class, which, using some specific technique, performs the time evolution.
 */
class Evolver {
public:
    virtual ~Evolver() = default;

    /**
     * @brief Prepares the evolution for a specific state, starting from t = 0 to @a maxTime.
     * @param initialState state to evolve
     * @param maxTime target time to evolve to
     * @param numSteps steps to divide the target time into. t = 0 is counted as the first step
     */
    virtual void prepareFor(const arma::cx_vec &initialState, double maxTime, std::size_t numSteps) = 0;

    /**
     * @brief Performs the evolution to the next step, increasing time as calculated in prepareFor().
     * @details The evolution is possible to at least one step past @a maxTime from prepareFor().
     */
    virtual void evolve() = 0;

    [[nodiscard]] virtual const arma::cx_vec &getCurrentState() const = 0;
};

#endif //MBL_ED_EVOLVER_H
