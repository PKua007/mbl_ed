//
// Created by Piotr Kubala on 22/01/2020.
//

#include "utils/Assertions.h"
#include "Eigensystem.h"

#include <utility>

Eigensystem::Eigensystem(const arma::vec &eigenvalues, const arma::mat &eigenvectors) {
    std::size_t size = eigenvalues.size();
    Expects(eigenvectors.n_cols == size);
    Expects(eigenvectors.n_rows == size);

    this->eigenenergies.reserve(size);
    this->eigenstates.resize(size);

    std::copy(eigenvalues.begin(), eigenvalues.end(), std::back_inserter(this->eigenenergies));
    for (std::size_t i = 0; i < size; i++)
        std::copy(eigenvectors.begin_col(i), eigenvectors.end_col(i), std::back_inserter(this->eigenstates[i]));
}

const std::vector<double> &Eigensystem::getEigenenergies() const {
    return this->eigenenergies;
}

const std::vector<std::vector<double>> &Eigensystem::getEigenstates() const {
    return this->eigenstates;
}

void Eigensystem::addEntry(double eigenenergy, const std::vector<double> &eigenvector) {
    if (!this->eigenstates.empty()) {
        const auto &firstState = this->eigenstates.front();
        Expects(eigenvector.size() == firstState.size());

        if (!firstState.empty() && this->isComplete())
            throw std::runtime_error("The number of eigenentries is already exhausted");
    }

    this->eigenenergies.push_back(eigenenergy);
    this->eigenstates.push_back(eigenvector);
}

bool Eigensystem::isComplete() const {
    if (this->eigenstates.empty())
        return true;

    const auto &firstState = this->eigenstates.front();
    if (firstState.empty())
        return true;

    return firstState.size() == this->eigenenergies.size();
}

std::size_t Eigensystem::size() const {
    return this->eigenenergies.size();
}

bool Eigensystem::empty() const {
    return this->eigenenergies.empty();
}

bool Eigensystem::hasEigenvectors() const {
    if (this->empty())
        return false;
    return !this->eigenstates.front().empty();
}

Eigensystem::Eigensystem(std::vector<double> eigenenergies) : eigenenergies{std::move(eigenenergies)} {
    this->eigenstates.resize(this->eigenenergies.size());
}

bool operator==(const Eigensystem &lhs, const Eigensystem &rhs) {
    return lhs.eigenenergies == rhs.eigenenergies &&
           lhs.eigenstates == rhs.eigenstates;
}

bool operator!=(const Eigensystem &lhs, const Eigensystem &rhs) {
    return !(rhs == lhs);
}
