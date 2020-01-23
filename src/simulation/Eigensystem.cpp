//
// Created by Piotr Kubala on 22/01/2020.
//

#include <utility>
#include <iterator>

#include <ZipIterator.hpp>

#include "utils/Assertions.h"
#include "Eigensystem.h"

Eigensystem::Eigensystem(arma::vec eigenvalues, arma::mat eigenvectors)
        : eigenenergies{std::move(eigenvalues)}, eigenstates{std::move(eigenvectors)}, hasEigenvectors_{true}
{
    std::size_t size = this->eigenenergies.size();
    Expects(this->eigenstates.n_cols == size);
    Expects(this->eigenstates.n_rows == size);
    for (std::size_t i{}; i < size; i++)
        Expects(arma::any(this->eigenstates.col(i)));

    if (size == 0)
        this->hasEigenvectors_ = false;
    this->sortAndNormalize();
}

std::size_t Eigensystem::size() const {
    return this->eigenenergies.size();
}

bool Eigensystem::empty() const {
    return this->eigenenergies.empty();
}

bool Eigensystem::hasEigenvectors() const {
    return this->hasEigenvectors_;
}

const arma::vec &Eigensystem::getEigenenergies() const {
    return this->eigenenergies;
}

const arma::mat &Eigensystem::getEigenstates() const {
    return this->eigenstates;
}

arma::vec Eigensystem::getEigenstate(std::size_t i) const {
    Expects(i < this->size());
    return this->eigenstates.col(i);
}

arma::vec Eigensystem::getNormalizedEigenenergies() const {
    if (this->empty())
        return {};
    else if (this->size() == 1)
        return {1};
    else if (arma::all(this->eigenenergies == this->eigenenergies.front()))
        throw std::runtime_error("All eigenvalues equal, cannot normalize.");

    double low = this->eigenenergies.front();
    double high = this->eigenenergies.back();
    return (this->eigenenergies - low) / (high - low);
}

void Eigensystem::store(std::ostream &eigenenergiesOut) const {
    if (!this->eigenenergies.save(eigenenergiesOut))
        throw std::runtime_error("Store procedure failed");
}

void Eigensystem::restore(std::istream &in) {
    arma::vec newEigenenergies;
    if (!newEigenenergies.load(in))
        throw std::runtime_error("Restore procedure failed");
    *this = Eigensystem(newEigenenergies);
}

Eigensystem::Eigensystem(const arma::vec &eigenvalues)
        : eigenenergies{eigenvalues}, eigenstates(eigenvalues.size(), eigenvalues.size(), arma::fill::zeros),
          hasEigenvectors_{false}
{
    this->sortAndNormalize();
}

bool operator==(const Eigensystem &lhs, const Eigensystem &rhs) {
    return arma::approx_equal(lhs.eigenenergies, rhs.eigenenergies, "absdiff", 1e-12) &&
           arma::approx_equal(lhs.eigenstates, rhs.eigenstates, "absdiff", 1e-12);
}

bool operator!=(const Eigensystem &lhs, const Eigensystem &rhs) {
    return !(rhs == lhs);
}

void Eigensystem::sortAndNormalize() {
    auto indices = arma::regspace<arma::ivec>(0, this->size() - 1);
    auto zipped = Zip(this->eigenenergies, indices);
    std::sort(zipped.begin(), zipped.end());

    if (!this->hasEigenvectors_)
        return;

    arma::mat newEigenstates(this->size(), this->size());
    for (std::size_t i{}; i < this->size(); i++)
        newEigenstates.col(i) = arma::normalise(this->eigenstates.col(indices[i]));
    this->eigenstates = newEigenstates;
}

std::ostream &operator<<(std::ostream &out, const Eigensystem &eigensystem) {
    out << eigensystem.eigenenergies.t() << std::endl << eigensystem.eigenstates;
    return out;
}
