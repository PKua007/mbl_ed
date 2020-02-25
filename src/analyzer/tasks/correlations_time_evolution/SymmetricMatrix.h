//
// Created by Piotr Kubala on 19/02/2020.
//

#ifndef MBL_ED_SYMMETRICMATRIX_H
#define MBL_ED_SYMMETRICMATRIX_H

#include <vector>
#include <armadillo>

#include "utils/Assertions.h"

template<typename T>
class SymmetricMatrix {
private:
    std::size_t size_{};
    std::vector<T> elements;

    [[nodiscard]] std::size_t idx(std::size_t i, std::size_t j) const {
        Expects(i < this->size_);
        Expects(j < this->size_);

        std::size_t smaller = std::min(i, j);
        std::size_t bigger = std::max(i, j);

        return bigger*(bigger + 1)/2 + smaller;
    }

public:
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;

    SymmetricMatrix() = default;
    explicit SymmetricMatrix(std::size_t size_) : size_{size_}, elements(size_*(size_ + 1)/2) { }

    [[nodiscard]] std::size_t size() const { return this->size_; }
    [[nodiscard]] T &operator()(std::size_t i, std::size_t j) { return this->elements[this->idx(i, j)]; }
    [[nodiscard]] const T &operator()(std::size_t i, std::size_t j) const { return this->elements[this->idx(i, j)]; }

    template<typename T_ = T>
    [[nodiscard]] typename std::enable_if_t<std::is_convertible_v<T_, double>, arma::mat> toArma() const {
        arma::mat arma_mat(this->size_, this->size_);
        for (std::size_t i = 0; i < this->size_; i++)
            for (std::size_t j = 0; j < this->size_; j++)
                arma_mat(i, j) = this->operator()(i, j);
        return arma_mat;
    }

    template<typename T_ = T>
    friend typename std::enable_if_t<std::is_convertible_v<T_, double>, std::ostream &>
    operator<<(std::ostream &out, const SymmetricMatrix &matrix) {
        return out << matrix.toArma();
    }

    [[nodiscard]] iterator begin() { return this->elements.begin(); }
    [[nodiscard]] iterator end() { return this->elements.end(); }
    [[nodiscard]] const_iterator begin() const { return this->elements.begin(); }
    [[nodiscard]] const_iterator end() const { return this->elements.end(); }
};

#endif //MBL_ED_SYMMETRICMATRIX_H
