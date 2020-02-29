//
// Created by Piotr Kubala on 19/02/2020.
//

#ifndef MBL_ED_SYMMETRICMATRIX_H
#define MBL_ED_SYMMETRICMATRIX_H

#include <vector>
#include <armadillo>

#include "utils/Assertions.h"

/**
 * @brief A class representing 2D symmetric square matrix, with elements of any type. The elements lying symmetrically
 * wrt. diagonal are actually the same objects in memory, which saves the space.
 * @tparam T type of the elements stored in the matrix. It may not be numeric.
 */
template<typename T>
class SymmetricMatrix {
private:
    std::size_t size_{};
    std::vector<T> elements;

    /**
     * @brief Converts the @a row and @a col to the index in internal storage taking into account symmetricity.
     */
    [[nodiscard]] std::size_t idx(std::size_t row, std::size_t col) const {
        Expects(row < this->size_);
        Expects(col < this->size_);

        std::size_t smaller = std::min(row, col);
        std::size_t bigger = std::max(row, col);

        return bigger*(bigger + 1)/2 + smaller;
    }

public:
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    SymmetricMatrix() = default;
    explicit SymmetricMatrix(std::size_t size_) : size_{size_}, elements(size_*(size_ + 1)/2) { }

    /**
     * @brief Returns the linear size of the matrix, ie. number of rows or columns.
     */
    [[nodiscard]] std::size_t size() const { return this->size_; }

    /**
     * @brief Mutable access to the elements. Elements are counted from 0. As it is the symmetric matrix, @a row and
     * @a col can always be swapped and this leads to the same object. The modifications also effectively modify both
     * permutations.
     */
    [[nodiscard]] T &operator()(std::size_t row, std::size_t col) { return this->elements[this->idx(row, col)]; }

    /**
     * @brief Immutable access to the elements. Elements are counted from 0. As it is the symmetric matrix, @a row and
     * @a col can always be swapped and this leads to the same object.
     */
    [[nodiscard]] const T &operator()(std::size_t row, std::size_t col) const {
        return this->elements[this->idx(row, col)];
    }

    /**
     * @brief A cast to arma dense double matrix, enabled if @a T is convertible to double.
     */
    template<typename T_ = T>
    [[nodiscard]] typename std::enable_if_t<std::is_convertible_v<T_, double>, arma::mat> toArma() const {
        arma::mat arma_mat(this->size_, this->size_);
        for (std::size_t i = 0; i < this->size_; i++)
            for (std::size_t j = 0; j < this->size_; j++)
                arma_mat(i, j) = this->operator()(i, j);
        return arma_mat;
    }

    /**
     * @brief Stream insertion operator, enabled if @a T is convertible to double (as it uses arma to print the matrix).
     */
    template<typename T_ = T>
    friend typename std::enable_if_t<std::is_convertible_v<T_, double>, std::ostream &>
    operator<<(std::ostream &out, const SymmetricMatrix &matrix) {
        return out << matrix.toArma();
    }

    /**
     * @brief Begin iterator. Iterators effectively iterate only the upper-triangular (or lower-triangular) part of the
     * matrix as a reflection of the symmetricity.
     */
    [[nodiscard]] iterator begin() { return this->elements.begin(); }
    [[nodiscard]] iterator end() { return this->elements.end(); }
    [[nodiscard]] const_iterator begin() const { return this->elements.begin(); }
    [[nodiscard]] const_iterator end() const { return this->elements.end(); }
};

#endif //MBL_ED_SYMMETRICMATRIX_H
