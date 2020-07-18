//
// Created by Piotr Kubala on 18/07/2020.
//

#ifndef MBL_ED_FOCKVECTOR_H
#define MBL_ED_FOCKVECTOR_H

#include <vector>

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

#include <vector>
#include <map>

#endif //MBL_ED_FOCKVECTOR_H
