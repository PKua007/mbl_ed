//
// Created by Piotr Kubala on 03/10/2020.
//

#ifndef MBL_ED_ARMAAPPROXEQUALTROMPELOEILMATCHER_H
#define MBL_ED_ARMAAPPROXEQUALTROMPELOEILMATCHER_H

#include <catch2/trompeloeil.hpp>
#include <armadillo>

inline auto arma_eq(const arma::cx_vec &vec) {
    return trompeloeil::make_matcher<arma::cx_vec>(
            [](const arma::cx_vec &value, const arma::cx_vec &expected) {
                return arma::approx_equal(value, expected, "absdiff", 1e-8);
            },
            [](std::ostream &os, const arma::cx_vec &expected) {
                os << " arma::approx_equal " << expected;
            },
            vec
    );
}

#endif //MBL_ED_ARMAAPPROXEQUALTROMPELOEILMATCHER_H
