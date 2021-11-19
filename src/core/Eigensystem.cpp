//
// Created by Piotr Kubala on 22/01/2020.
//

#include <utility>
#include <iterator>

#include <ZipIterator.hpp>

#include "utils/Assertions.h"
#include "Eigensystem.h"

Eigensystem::Eigensystem(arma::vec eigenvalues, arma::mat eigenstates, std::shared_ptr<const FockBasis> fockBasis)
        : eigenenergies{std::move(eigenvalues)}, eigenstates{std::move(eigenstates)}, hasEigenvectors_{true},
          fockBasis{std::move(fockBasis)}
{
    std::size_t size = this->eigenenergies.size();
    if (this->fockBasis != nullptr)
        Expects(this->fockBasis->size() == size);
    Expects(this->eigenstates.n_cols == size);
    Expects(this->eigenstates.n_rows == size);
    for (std::size_t i{}; i < size; i++)
        Expects(arma::any(this->eigenstates.col(i)));

    if (size == 0)
        this->hasEigenvectors_ = false;
    this->sortEigenenergiesAndNormalizeEigenstates();
}

Eigensystem::Eigensystem(arma::vec eigenvalues, std::shared_ptr<const FockBasis> fockBasis)
        : eigenenergies{std::move(eigenvalues)}, hasEigenvectors_{false}, fockBasis{std::move(fockBasis)}
{
    if (this->fockBasis != nullptr)
        Expects(this->fockBasis->size() == this->eigenenergies.size());
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
        throw std::runtime_error("Eigenenergies store procedure failed");
}

void Eigensystem::store(std::ostream &eigenenergiesOut, std::ostream &eigenstatesOut) const {
    if (!this->eigenenergies.save(eigenenergiesOut))
        throw std::runtime_error("Eigenenergies store procedure failed");
    if (!this->eigenstates.save(eigenstatesOut))
        throw std::runtime_error("Eigenstates store procedure failed");
}

void Eigensystem::restore(std::istream &eigenenergiesIn, std::shared_ptr<const FockBasis> newFockBasis) {
    arma::vec newEigenenergies;
    if (!newEigenenergies.load(eigenenergiesIn))
        throw std::runtime_error("Eigenenergies restore procedure failed");
    *this = Eigensystem(newEigenenergies, std::move(newFockBasis));
}

void Eigensystem::restore(std::istream &eigenenergiesIn, std::istream &eigenstatesIn,
                          std::shared_ptr<const FockBasis> newFockBasis)
{
    arma::vec newEigenenergies;
    arma::mat newEigenstates;
    if (!newEigenenergies.load(eigenenergiesIn))
        throw std::runtime_error("Eigenenergies restore procedure failed");
    if (!newEigenstates.load(eigenstatesIn))
        throw std::runtime_error("Eigenstates restore procedure failed");
    *this = Eigensystem(newEigenenergies, newEigenstates, std::move(newFockBasis));
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
    // We let the endpoints of the range to be outside [0, 1], because it's hard to catch 0 or 1 otherwise due to
    // machine precision - we no longer use Expects(epsilon - delta/2 >= 0 && epsilon + delta/2 <= 1)

    auto normalizedEnergies = this->getNormalizedEigenenergies();
    double relativeFrom = epsilon - delta/2;
    double relativeTo = epsilon + delta/2;

    auto fromIt = std::lower_bound(normalizedEnergies.begin(), normalizedEnergies.end(), relativeFrom);
    auto toIt = std::lower_bound(normalizedEnergies.begin(), normalizedEnergies.end(), relativeTo);

    std::vector<std::size_t> indices(toIt - fromIt);
    std::iota(indices.begin(), indices.end(), fromIt - normalizedEnergies.begin());
    return indices;
}

std::vector<std::size_t> Eigensystem::getIndicesOfNumberOfNormalizedEnergies(double epsilon,
                                                                             std::size_t numEnergies) const
{
    Expects(epsilon > 0 && epsilon < 1);
    Expects(numEnergies > 0);
    Expects(numEnergies <= this->eigenenergies.size());

    auto normalizedEnergies = this->getNormalizedEigenenergies();
    std::vector<double> distances(normalizedEnergies.size());
    std::transform(normalizedEnergies.begin(), normalizedEnergies.end(), distances.begin(),
                   [epsilon](double e) { return std::abs(e - epsilon); });

    std::vector<std::size_t> indices(this->eigenenergies.size());
    std::iota(indices.begin(), indices.end(), 0);

    auto distZipped = Zip(distances, indices);
    std::sort(distZipped.begin(), distZipped.end());

    indices.resize(numEnergies);
    std::sort(indices.begin(), indices.end());
    return indices;
}

const FockBasis &Eigensystem::getFockBasis() const {
    return *(this->fockBasis);
}

bool Eigensystem::hasFockBasis() const {
    return this->fockBasis != nullptr;
}

bool Eigensystem::isOrthonormal() const {
    if (!this->hasEigenvectors_)
        throw std::runtime_error("Eigensystem does not contain eigenvectors");

    // we have up to this->fockBasis->size() additions
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
