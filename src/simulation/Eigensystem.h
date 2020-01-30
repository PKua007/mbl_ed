//
// Created by Piotr Kubala on 22/01/2020.
//

#ifndef MBL_ED_EIGENSYSTEM_H
#define MBL_ED_EIGENSYSTEM_H

#include <vector>
#include <armadillo>

class Eigensystem {
private:
    arma::vec eigenenergies;
    arma::mat eigenstates;
    bool hasEigenvectors_{};

    void sortAndNormalize();
    void sort();

public:
    Eigensystem() = default;
    Eigensystem(arma::vec eigenvalues, arma::mat eigenvectors);
    explicit Eigensystem(arma::vec eigenvalues);

    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool hasEigenvectors() const;
    [[nodiscard]] const arma::vec &getEigenenergies() const;
    [[nodiscard]] const arma::mat &getEigenstates() const;
    [[nodiscard]] arma::vec getEigenstate(std::size_t i) const;
    [[nodiscard]] arma::vec getNormalizedEigenenergies() const;
    [[nodiscard]] std::vector<std::size_t> getIndicesOfNormalizedEnergiesInBand(double epsilon, double delta) const;
    void store(std::ostream &eigenenergiesOut) const;
    void restore(std::istream &in);

    friend bool operator==(const Eigensystem &lhs, const Eigensystem &rhs);
    friend bool operator!=(const Eigensystem &lhs, const Eigensystem &rhs);
    friend std::ostream &operator<<(std::ostream &out, const Eigensystem &eigensystem);
};

#endif //MBL_ED_EIGENSYSTEM_H
