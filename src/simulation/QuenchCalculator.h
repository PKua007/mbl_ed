//
// Created by Piotr Kubala on 18/08/2020.
//

#ifndef MBL_ED_QUENCHCALCULATOR_H
#define MBL_ED_QUENCHCALCULATOR_H

#include <vector>
#include <armadillo>


class QuenchCalculator {
private:
    static constexpr std::size_t MIN_EIGVALS = 6;

    std::vector<double> quenchEpsilons;
    std::vector<double> quenchEpsilonVariances;

    arma::vec lastQuenchedState;

public:
    void addQuench(const arma::sp_mat &initialHamiltonian, const arma::sp_mat &finalHamiltonian);

    [[nodiscard]] arma::vec getLastQuenchedState() const { return lastQuenchedState; }
    [[nodiscard]] double getLastQuenchEpsilon() const;
    [[nodiscard]] double getLastQuenchEpsilonQuantumUncertainty() const;

    [[nodiscard]] double getMeanEpsilon() const;
    [[nodiscard]] double getMeanEpsilonQuantumUncertainty() const;
    [[nodiscard]] double getEpsilonAveragingSampleError() const;
};


#endif //MBL_ED_QUENCHCALCULATOR_H
