//
// Created by pkua on 01.11.2019.
//

#ifndef MBL_ED_FOCKBASE_H
#define MBL_ED_FOCKBASE_H


#include <vector>
#include <map>

class FockVector {
private:
    std::vector<int> data;

public:
    using iterator = std::vector<int>::iterator;
    using const_iterator = std::vector<int>::const_iterator;

    FockVector() = default;
    explicit FockVector(std::size_t size, int init = 0) : data(size, init) { }
    FockVector(const std::initializer_list<int> &list) : data(list) { }

    [[nodiscard]] bool empty() const { return data.empty(); }
    [[nodiscard]] std::size_t size() const { return data.size(); }
    [[nodiscard]] int &operator[](std::size_t idx) { return data[idx]; }
    [[nodiscard]] int operator[](std::size_t idx) const { return data[idx]; }

    [[nodiscard]] int front() const { return data.front(); }
    [[nodiscard]] int back() const { return data.back(); }
    [[nodiscard]] int &front() { return data.front(); }
    [[nodiscard]] int &back() { return data.back(); }

    [[nodiscard]] iterator begin() { return data.begin(); }
    [[nodiscard]] iterator end() { return data.end(); }
    [[nodiscard]] const_iterator begin() const { return data.begin(); }
    [[nodiscard]] const_iterator end() const { return data.end(); }

    [[nodiscard]] friend bool operator==(const FockVector &fw1, const FockVector &fw2) { return fw1.data == fw2.data; }
    [[nodiscard]] friend bool operator!=(const FockVector &fw1, const FockVector &fw2) { return fw1.data == fw2.data; }
};

/**
 * @brief A class representing a base of product states of bosons/fermion trapped inside an optical lattice.
 * @details It uses hashing technique for a fast access to the elements.
 */
class FockBase {
public:
    using Vector = FockVector;
    using iterator = std::vector<Vector>::iterator;
    using const_iterator = std::vector<Vector>::const_iterator;

private:
    std::vector<Vector> theBase;
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


#endif //MBL_ED_FOCKBASE_H
