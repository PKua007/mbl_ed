//
// Created by Piotr Kubala on 21/02/2020.
//

#ifndef MBL_ED_VECTORAPPROXEQUALMATCHER_H
#define MBL_ED_VECTORAPPROXEQUALMATCHER_H

#include <catch2/catch.hpp>
#include <sstream>
#include <vector>
#include <iterator>

template<typename T>
class VectorApproxEqualMatcher : public Catch::MatcherBase<std::vector<T>> {
private:
    std::vector<T> expected;
    double epsilon;

public:
    VectorApproxEqualMatcher(std::vector<T> expected, double epsilon) : expected(std::move(expected)), epsilon(epsilon) { }

    bool match(const std::vector<T> &actual) const override {
        if (this->expected.size() != actual.size())
            return false;
        return std::equal(actual.begin(), actual.end(), this->expected.begin(),
                          [this](auto d1, auto d2) { return d1 == Approx(d2).epsilon(epsilon); });
    }

    [[nodiscard]] std::string describe() const override {
        std::ostringstream ss;
        ss << "is, within " << epsilon << " tolerance threshold, equal to" << std::endl;
        std::copy(this->expected.begin(), this->expected.end(), std::ostream_iterator<T>(ss, " "));
        return ss.str();
    }
};

template<typename T>
inline VectorApproxEqualMatcher<T> IsApproxEqual(const std::vector<T> &expected, double epsilon) {
    return VectorApproxEqualMatcher(expected, epsilon);
}

#endif //MBL_ED_VECTORAPPROXEQUALMATCHER_H
