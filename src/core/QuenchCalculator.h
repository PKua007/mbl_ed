//
// Created by Piotr Kubala on 18/08/2020.
//

#ifndef MBL_ED_QUENCHCALCULATOR_H
#define MBL_ED_QUENCHCALCULATOR_H

#include <vector>
#include <armadillo>

#include "simulation/Restorable.h"

/**
 * @brief A class performing quantum quenches from initial to final Hamiltonian.
 * @details A quench is when a ground state of initial Hamiltonian is considered w.r.t. final HamiltonianAfter.
 * The class calculates its normalized energy and quantum spread. The class can accept multiple quenches, from which
 * it calculates mean energy, its sample spread and mean quantum uncertainty.
 */
class QuenchCalculator : public Restorable {
private:
    static constexpr std::size_t MIN_EIGVALS = 6;

    std::vector<double> quenchEpsilons;
    std::vector<double> quenchEpsilonVariances;

    arma::vec lastQuenchedState;

public:
    /**
     * @brief Adds another quench to ensemble.
     */
    void addQuench(const arma::sp_mat &initialHamiltonian, const arma::sp_mat &finalHamiltonian);

    /**
     * @brief Returns the ground state of initial Hamiltonian from last addQuench invocation.
     */
    [[nodiscard]] const arma::vec &getLastQuenchedState() const { return lastQuenchedState; }

    /**
     * @brief Returns the normalized energy of ground state of initial Hamiltonian w.r.t. final Hamiltonian from last
     * addQuench invocation.
     */
    [[nodiscard]] double getLastQuenchEpsilon() const;

    /**
     * @brief Returns the quantum spread in normalized energy of ground state of initial Hamiltonian w.r.t. final
     * Hamiltonian from last addQuench invocation.
     */
    [[nodiscard]] double getLastQuenchEpsilonQuantumUncertainty() const;

    /**
     * @brief Returns mean normalized energy of a quenched state in the whole ensemble.
     */
    [[nodiscard]] double getMeanEpsilon() const;

    /**
     * @brief Returns mean quantum spread in normalized energy of a quenched state in the whole ensemble calculated
     * as a square root of mean quantum variances.
     */
    [[nodiscard]] double getMeanEpsilonQuantumUncertainty() const;

    /**
     * @brief Returns sample error when calculating mean normalized energy of the ensemble - it is the standard
     * deviation of the distribution of quench energies.
     */
    [[nodiscard]] double getEpsilonAveragingSampleError() const;

    void clear() override;
    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
};


#endif //MBL_ED_QUENCHCALCULATOR_H
