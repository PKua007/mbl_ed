//
// Created by pkua on 15.11.2019.
//

#ifndef MBL_ED_ARMAMATRIXAPPROXEQUALMATCHER_H
#define MBL_ED_ARMAMATRIXAPPROXEQUALMATCHER_H

#include <Catch2/catch.hpp>
#include <armadillo>
#include <ostream>
#include <utility>

class ArmaMatApproxEqualMatcher : public Catch::MatcherBase<arma::mat> {
private:
    arma::mat expected;
    double epsilon;

public:
    ArmaMatApproxEqualMatcher(arma::mat expected, double epsilon) : expected(std::move(expected)), epsilon(epsilon) { }

    bool match(const arma::mat &actual) const override {
        return arma::approx_equal(actual, expected, "absdiff", epsilon);
    }

    std::string describe() const override {
        std::ostringstream ss;
        ss << "is, within " << epsilon << " tolerance threshold, equal to" << std::endl << expected;
        return ss.str();
    }
};

inline ArmaMatApproxEqualMatcher IsApproxEqual(const arma::mat &expected, double epsilon) {
    return ArmaMatApproxEqualMatcher(expected, epsilon);
}

#endif //MBL_ED_ARMAMATRIXAPPROXEQUALMATCHER_H
