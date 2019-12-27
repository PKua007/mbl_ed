//
// Created by pkua on 01.11.2019.
//

#ifndef MBL_ED_FOCKBASE_H
#define MBL_ED_FOCKBASE_H


#include <vector>
#include <map>

class FockBase {
public:
    using Vector = std::vector<int>;
    using iterator = std::vector<Vector>::iterator;
    using const_iterator = std::vector<Vector>::const_iterator;

private:
    std::vector<Vector> theBase;
    std::map<double, std::size_t> indexMap;

    [[nodiscard]] double computeHash(const Vector &vector) const;

public:
    void add(Vector vector);
    [[nodiscard]] std::size_t size() const;
    Vector &operator[](std::size_t i);
    const Vector &operator[](std::size_t i) const;
    [[nodiscard]] std::optional<std::size_t> findIndex(const Vector &vector) const;
    iterator begin();
    iterator end();
    [[nodiscard]] const_iterator begin() const;
    [[nodiscard]] const_iterator end() const;
    [[nodiscard]] std::size_t getNumberOfSites() const;
    [[nodiscard]] std::size_t getNumberOfParticles() const;
};


#endif //MBL_ED_FOCKBASE_H
