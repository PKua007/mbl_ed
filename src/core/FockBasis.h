//
// Created by pkua on 01.11.2019.
//

#ifndef MBL_ED_FOCKBASIS_H
#define MBL_ED_FOCKBASIS_H


#include <vector>
#include <map>
#include <optional>

#include "FockVector.h"

/**
 * @brief A class representing a basis of product states of bosons/fermion trapped inside an optical lattice.
 * @details It uses hashing technique for a fast access to the elements.
 */
class FockBasis {
public:
    using Vector = FockVector;
    using iterator = std::vector<Vector>::iterator;
    using const_iterator = std::vector<Vector>::const_iterator;

private:
    std::vector<Vector> theBasis;
    std::map<double, std::size_t> indexMap;

    [[nodiscard]] double computeHash(const Vector &vector) const;

public:
    void add(Vector vector);
    [[nodiscard]] std::size_t size() const;

    /**
     * @brief Returns modifiable vector of index @a i. j-th elements of vector reperesents number of parcitles on site
     * j.
     */
    Vector &operator[](std::size_t i);

    /**
     * @brief Returns non-modifiable vector of index @a i. j-th elements of vector reperesents number of parcitles on
     * site j.
     */
    const Vector &operator[](std::size_t i) const;

    /**
     * @brief Returns the index of a given Vector or std::nullopt if it is not present.
     */
    [[nodiscard]] std::optional<std::size_t> findIndex(const Vector &vector) const;
    iterator begin();
    iterator end();
    [[nodiscard]] const_iterator begin() const;
    [[nodiscard]] const_iterator end() const;
    [[nodiscard]] std::size_t getNumberOfSites() const;
    [[nodiscard]] std::size_t getNumberOfParticles() const;
};


#endif //MBL_ED_FOCKBASIS_H
