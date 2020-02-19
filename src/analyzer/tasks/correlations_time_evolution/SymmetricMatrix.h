//
// Created by Piotr Kubala on 19/02/2020.
//

#ifndef MBL_ED_SYMMETRICMATRIX_H
#define MBL_ED_SYMMETRICMATRIX_H

#include <vector>

#include "utils/Assertions.h"

class SymmetricMatrix {
private:
    std::size_t size_{};
    std::vector<double> elements;

    [[nodiscard]] std::size_t idx(std::size_t i, std::size_t j) const {
        Expects(i < this->size_);
        Expects(j < this->size_);

        std::size_t smaller = std::min(i, j);
        std::size_t bigger = std::max(i, j);

        return bigger*(bigger + 1)/2 + smaller;
    }

public:
    SymmetricMatrix() = default;
    explicit SymmetricMatrix(std::size_t size_) : size_{size_}, elements(size_*(size_ + 1)/2) { }

    [[nodiscard]] std::size_t size() const { return this->size_; }
    [[nodiscard]] double &operator()(std::size_t i, std::size_t j) { return this->elements[this->idx(i, j)]; }
    [[nodiscard]] double operator()(std::size_t i, std::size_t j) const { return this->elements[this->idx(i, j)]; }
};


#endif //MBL_ED_SYMMETRICMATRIX_H
