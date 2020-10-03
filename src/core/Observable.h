//
// Created by Piotr Kubala on 18/09/2020.
//

#ifndef MBL_ED_OBSERVABLE_H
#define MBL_ED_OBSERVABLE_H

#include <vector>
#include <string>
#include <armadillo>

/**
 * @brief An interface representing an observable, or rather set of observables, which is calculated in some
 * implementation-dependent manner and its values can be accessed.
 */
class Observable {
public:
    /**
     * @brief Helper method calculating expected value of a diagonal @a observable or a @a state.
     */
    static double calculateExpectedValue(const arma::vec &observable, const arma::cx_vec &state) {
        arma::vec modulusSquare(state.size());
        for (std::size_t i = 0; i < state.size(); i++)
            modulusSquare[i] = state[i].real() * state[i].real() + state[i].imag() * state[i].imag();

        return arma::as_scalar(modulusSquare.t() * observable);
    }

    virtual ~Observable() = default;

    /**
     * @brief Returns the names of all observables in the same order as values.
     */
    [[nodiscard]] virtual std::vector<std::string> getHeader() const = 0;

    /**
     * @brief Returns the values of all observables in the same order as names.
     */
    [[nodiscard]] virtual std::vector<double> getValues() const = 0;
};


#endif //MBL_ED_OBSERVABLE_H
