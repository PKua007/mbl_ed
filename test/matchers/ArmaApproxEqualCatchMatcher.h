//
// Created by pkua on 15.11.2019.
//

#ifndef MBL_ED_ARMAAPPROXEQUALCATCHMATCHER_H
#define MBL_ED_ARMAAPPROXEQUALCATCHMATCHER_H

#include <catch2/catch.hpp>
#include <armadillo>
#include <ostream>
#include <utility>

template<typename T>
class ArmaApproxEqualMatcher : public Catch::MatcherBase<T> {
private:
    T expected;
    double epsilon;

public:
    ArmaApproxEqualMatcher(T expected, double epsilon) : expected(std::move(expected)), epsilon(epsilon) { }

    bool match(const T &actual) const override {
        return arma::approx_equal(actual, expected, "absdiff", epsilon);
    }

    [[nodiscard]] std::string describe() const override {
        std::ostringstream ss;
        ss << "is, within " << epsilon << " tolerance threshold, equal to" << std::endl << expected;
        return ss.str();
    }
};

template<typename T>
inline ArmaApproxEqualMatcher<T> IsApproxEqual(const T &expected, double epsilon) {
    return ArmaApproxEqualMatcher(expected, epsilon);
}

#endif //MBL_ED_ARMAAPPROXEQUALCATCHMATCHER_H
