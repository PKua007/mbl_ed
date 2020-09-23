//
// Created by Piotr Kubala on 18/09/2020.
//

#ifndef MBL_ED_OBSERVABLE_H
#define MBL_ED_OBSERVABLE_H

#include <vector>
#include <string>
#include <armadillo>

class Observable {
public:
    static double calculateExpectedValue(const arma::vec &observable, const arma::cx_vec &state) {
        arma::vec modulusSquare(state.size());
        for (std::size_t i = 0; i < state.size(); i++)
            modulusSquare[i] = state[i].real() * state[i].real() + state[i].imag() * state[i].imag();

        return arma::as_scalar(modulusSquare.t() * observable);
    }

    virtual ~Observable() = default;

    [[nodiscard]] virtual std::vector<std::string> getHeader() const = 0;
    [[nodiscard]] virtual std::vector<double> getValues() const = 0;
};


#endif //MBL_ED_OBSERVABLE_H
