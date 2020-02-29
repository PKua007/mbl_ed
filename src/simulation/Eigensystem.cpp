//
// Created by Piotr Kubala on 22/01/2020.
//

#include <utility>
#include <iterator>

#include <ZipIterator.hpp>

#include "utils/Assertions.h"
#include "Eigensystem.h"

Eigensystem::Eigensystem(arma::vec eigenvalues, arma::mat eigenvectors, std::shared_ptr<const FockBase> fockBase)
        : eigenenergies{std::move(eigenvalues)}, eigenstates{std::move(eigenvectors)}, hasEigenvectors_{true},
          fockBase{std::move(fockBase)}
{
    std::size_t size = this->eigenenergies.size();
    if (this->fockBase != nullptr)
        Expects(this->fockBase->size() == size);
    Expects(this->eigenstates.n_cols == size);
    Expects(this->eigenstates.n_rows == size);
    for (std::size_t i{}; i < size; i++)
        Expects(arma::any(this->eigenstates.col(i)));

    if (size == 0)
        this->hasEigenvectors_ = false;
    this->sortEigenenergiesAndNormalizeEigenstates();
}

Eigensystem::Eigensystem(arma::vec eigenvalues, std::shared_ptr<const FockBase> fockBase)
        : eigenenergies{std::move(eigenvalues)}, hasEigenvectors_{false}, fockBase{std::move(fockBase)}
{
    if (this->fockBase != nullptr)
        Expects(this->fockBase->size() == this->eigenenergies.size());
    this->sortEigenenergies();
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
    if (!this->hasEigenvectors_)
        throw std::runtime_error("Eigensystem does not contain eigenvectors");
    return this->eigenstates;
}

arma::vec Eigensystem::getEigenstate(std::size_t i) const {
    if (!this->hasEigenvectors_)
        throw std::runtime_error("Eigensystem does not contain eigenvectors");
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

void Eigensystem::restore(std::istream &in, std::shared_ptr<const FockBase> newFockBase) {
    arma::vec newEigenenergies;
    if (!newEigenenergies.load(in))
        throw std::runtime_error("Restore procedure failed");
    *this = Eigensystem(newEigenenergies, std::move(newFockBase));
}

bool operator==(const Eigensystem &lhs, const Eigensystem &rhs) {
    return arma::approx_equal(lhs.eigenenergies, rhs.eigenenergies, "absdiff", 1e-12) &&
           arma::approx_equal(lhs.eigenstates, rhs.eigenstates, "absdiff", 1e-12);
}

bool operator!=(const Eigensystem &lhs, const Eigensystem &rhs) {
    return !(rhs == lhs);
}

void Eigensystem::sortEigenenergiesAndNormalizeEigenstates() {
    auto indices = arma::regspace<arma::ivec>(0, this->size() - 1);
    auto zipped = Zip(this->eigenenergies, indices);
    std::sort(zipped.begin(), zipped.end());

    arma::mat newEigenstates(this->size(), this->size());
    for (std::size_t i{}; i < this->size(); i++)
        newEigenstates.col(i) = arma::normalise(this->eigenstates.col(indices[i]));
    this->eigenstates = newEigenstates;
}

void Eigensystem::sortEigenenergies() {
    std::sort(this->eigenenergies.begin(), this->eigenenergies.end());
}

std::ostream &operator<<(std::ostream &out, const Eigensystem &eigensystem) {
    out << eigensystem.eigenenergies.t() << std::endl << eigensystem.eigenstates;
    return out;
}

std::vector<std::size_t> Eigensystem::getIndicesOfNormalizedEnergiesInBand(double epsilon, double delta) const {
    Expects(epsilon > 0 && epsilon < 1);
    Expects(delta > 0);
    Expects(epsilon - delta/2 >= 0 && epsilon + delta/2 <= 1);

    auto normalizedEnergies = this->getNormalizedEigenenergies();
    double relativeFrom = epsilon - delta/2;
    double relativeTo = epsilon + delta/2;

    auto fromIt = std::lower_bound(normalizedEnergies.begin(), normalizedEnergies.end(), relativeFrom);
    auto toIt = std::lower_bound(normalizedEnergies.begin(), normalizedEnergies.end(), relativeTo);
    Assert(fromIt - normalizedEnergies.begin() >= 1);
    Assert(normalizedEnergies.end() - toIt >= 1);
    Assert(toIt - fromIt > 0);

    std::vector<std::size_t> indices(toIt - fromIt);
    std::iota(indices.begin(), indices.end(), fromIt - normalizedEnergies.begin());
    return indices;
}

const FockBase &Eigensystem::getFockBase() const {
    return *(this->fockBase);
}

bool Eigensystem::hasFockBase() const {
    return this->fockBase != nullptr;
}

bool Eigensystem::isOrthonormal() const {
    if (!this->hasEigenvectors_)
        throw std::runtime_error("Eigensystem does not contain eigenvectors");

    // we have up to this->fockBase->size() additions
    double epsilon = std::numeric_limits<double>::epsilon() * this->size();

    for (std::size_t i{}; i < this->size(); i++) {
        for (std::size_t j = i + 1; j < this->size(); j++) {
            double product = arma::dot(this->getEigenstate(i), this->getEigenstate(j));
            if (std::abs(product) > epsilon)
                return false;
        }
    }
    return true;
}
