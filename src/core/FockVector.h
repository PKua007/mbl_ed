//
// Created by Piotr Kubala on 18/07/2020.
//

#ifndef MBL_ED_FOCKVECTOR_H
#define MBL_ED_FOCKVECTOR_H

#include <vector>
#include <string>
#include <stdexcept>

struct FockVectorParseException : public std::runtime_error {
    explicit FockVectorParseException(const std::string &what) : std::runtime_error(what) { }
};

class FockVector {
private:
    std::vector<int> data;

public:
    using iterator = std::vector<int>::iterator;
    using const_iterator = std::vector<int>::const_iterator;

    FockVector() = default;
    explicit FockVector(std::size_t size, int init = 0) : data(size, init) { }
    FockVector(const std::initializer_list<int> &list) : data(list) { }
    explicit FockVector(const std::string &occupationRepresentation);
    FockVector(std::size_t sites, const std::string &tag);

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

    /**
     * @brief Concatenated @a fw1 and @a fw2 into one, larger vector.
     */
    [[nodiscard]] friend FockVector operator+(const FockVector &fw1, const FockVector &fw2);

    [[nodiscard]] friend bool operator==(const FockVector &fw1, const FockVector &fw2) { return fw1.data == fw2.data; }
    [[nodiscard]] friend bool operator!=(const FockVector &fw1, const FockVector &fw2) { return fw1.data == fw2.data; }

    friend std::ostream &operator<<(std::ostream &out, const FockVector &fw);
};

#endif //MBL_ED_FOCKVECTOR_H
