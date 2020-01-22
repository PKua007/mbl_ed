//
// Created by Piotr Kubala on 22/01/2020.
//

#ifndef MBL_ED_EIGENSYSTEM_H
#define MBL_ED_EIGENSYSTEM_H

#include <vector>
#include <armadillo>

class Eigensystem {
private:
    std::vector<double> eigenenergies;
    std::vector<std::vector<double>> eigenstates;

public:
    Eigensystem() = default;
    Eigensystem(std::vector<double> eigenenergies);
    Eigensystem(const arma::vec &eigenvalues, const arma::mat &eigenvectors);

    void addEntry(double eigenenergy, const std::vector<double> &eigenvector);
    [[nodiscard]] bool empty() const;
    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] const std::vector<double> &getEigenenergies() const;
    [[nodiscard]] const std::vector<std::vector<double>> &getEigenstates() const;
    [[nodiscard]] bool isComplete() const;
    [[nodiscard]] bool hasEigenvectors() const;

    friend bool operator==(const Eigensystem &lhs, const Eigensystem &rhs);
    friend bool operator!=(const Eigensystem &lhs, const Eigensystem &rhs);
};

#endif //MBL_ED_EIGENSYSTEM_H
