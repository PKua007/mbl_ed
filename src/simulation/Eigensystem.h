//
// Created by Piotr Kubala on 22/01/2020.
//

#ifndef MBL_ED_EIGENSYSTEM_H
#define MBL_ED_EIGENSYSTEM_H

#include <vector>
#include <armadillo>

#include "FockBase.h"

/**
 * @brief A system of eigenenergies in ascending order with corresponding eigenvectors normalized to a unity.
 *
 * The system optionaly may not contain eigenvectors.
 */
class Eigensystem {
private:
    arma::vec eigenenergies;
    arma::mat eigenstates;
    bool hasEigenvectors_{};
    std::shared_ptr<const FockBase> fockBase;

    void sortEigenenergiesAndNormalizeEigenstates();
    void sortEigenenergies();

public:
    Eigensystem() = default;

    /**
     * @brief Constructs a system: eigenvalues are entries in @a eigenvalues vector and eigenvectors are corresponding
     * columns in @a eigenvectors matrix.
     *
     * Eigenvalues are sorted in ascending order and eigenvectors are normalized to unity. @a fockBase can be omitted,
     * but when it's not, the size must match the number of eigenvalues.
     */
    Eigensystem(arma::vec eigenvalues, arma::mat eigenvectors, std::shared_ptr<const FockBase> fockBase = nullptr);

    /**
     * @brief Constructs a system without eigenvectors: eigenvalues are entries in @a eigengenvalues.
     *
     * Eigenvalues are sorted in ascending order.
     */
    explicit Eigensystem(arma::vec eigenvalues, std::shared_ptr<const FockBase> fockBase = nullptr);

    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool hasEigenvectors() const;
    [[nodiscard]] const arma::vec &getEigenenergies() const;
    [[nodiscard]] const arma::mat &getEigenstates() const;
    [[nodiscard]] bool hasFockBase() const;
    [[nodiscard]] const FockBase &getFockBase() const;
    [[nodiscard]] bool isOrthonormal() const;

    /**
     * @brief Returns eigenstate as column vector for eigenenergy of index @a i (in ascending order)
     */
    [[nodiscard]] arma::vec getEigenstate(std::size_t i) const;

    /**
     * @brief Returns eigenenrgies in the ascending order, but linearly normalized to be from [0, 1]
     */
    [[nodiscard]] arma::vec getNormalizedEigenenergies() const;

    /**
     * @brief Returns indices in vector from getNormalizedEigenenergies() corresponding to energies from a band
     * specified  by @a epsilon and @a delta.
     * @param epsilon the middle energy from the band (from [0, 1])
     * @param delta the width of the band
     */
    [[nodiscard]] std::vector<std::size_t> getIndicesOfNormalizedEnergiesInBand(double epsilon, double delta) const;
    void store(std::ostream &eigenenergiesOut) const;
    void restore(std::istream &in, std::shared_ptr<const FockBase> fockBase = nullptr);

    friend bool operator==(const Eigensystem &lhs, const Eigensystem &rhs);
    friend bool operator!=(const Eigensystem &lhs, const Eigensystem &rhs);
    friend std::ostream &operator<<(std::ostream &out, const Eigensystem &eigensystem);
};

#endif //MBL_ED_EIGENSYSTEM_H
