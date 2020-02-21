//
// Created by Piotr Kubala on 20/02/2020.
//

#include <catch2/catch.hpp>

#include "matchers/ArmaApproxEqualMatcher.h"

#include "analyzer/tasks/correlations_time_evolution/SymmetricMatrix.h"

TEST_CASE("SymmetricMatrix: default initialization") {
    SECTION("empty") {
        SymmetricMatrix matrix{};

        Expects(matrix.size() == 0);
    }

    SECTION("non-empty") {
        SymmetricMatrix matrix(3);

        Expects(matrix.size() == 3);
        Expects(matrix(0, 0) == 0); Expects(matrix(0, 1) == 0); Expects(matrix(0, 2) == 0);
        Expects(matrix(1, 0) == 0); Expects(matrix(1, 1) == 0); Expects(matrix(1, 2) == 0);
        Expects(matrix(2, 0) == 0); Expects(matrix(2, 1) == 0); Expects(matrix(2, 2) == 0);
    }
}

TEST_CASE("SymmetricMatrix: read write") {
    SECTION("indices: smaller, bigger") {
        SymmetricMatrix matrix(3);

        matrix(0, 0) = 1;
        matrix(0, 1) = 2;
        matrix(0, 2) = 3;
        matrix(1, 1) = 4;
        matrix(1, 2) = 5;
        matrix(2, 2) = 6;

        Expects(matrix(0, 0) == 1); Expects(matrix(0, 1) == 2); Expects(matrix(0, 2) == 3);
        Expects(matrix(1, 0) == 2); Expects(matrix(1, 1) == 4); Expects(matrix(1, 2) == 5);
        Expects(matrix(2, 0) == 3); Expects(matrix(2, 1) == 5); Expects(matrix(2, 2) == 6);
    }

    SECTION("indices: bigger, smaller") {
        SymmetricMatrix matrix(3);

        matrix(0, 0) = 1;
        matrix(1, 0) = 2;
        matrix(2, 0) = 3;
        matrix(1, 1) = 4;
        matrix(2, 1) = 5;
        matrix(2, 2) = 6;

        Expects(matrix(0, 0) == 1); Expects(matrix(0, 1) == 2); Expects(matrix(0, 2) == 3);
        Expects(matrix(1, 0) == 2); Expects(matrix(1, 1) == 4); Expects(matrix(1, 2) == 5);
        Expects(matrix(2, 0) == 3); Expects(matrix(2, 1) == 5); Expects(matrix(2, 2) == 6);
    }
}

TEST_CASE("SymmetricMatrix: bounds check") {
    SymmetricMatrix matrix(3);

    REQUIRE_THROWS(matrix(0, 3));
    REQUIRE_THROWS(matrix(3, 0));
    REQUIRE_THROWS(matrix(4, 1));
    REQUIRE_THROWS(matrix(1, 4));
}

TEST_CASE("SymmetricMatrix: toArma") {
    SymmetricMatrix matrix(3);

    matrix(0, 0) = 1;
    matrix(1, 0) = 2;
    matrix(2, 0) = 3;
    matrix(1, 1) = 4;
    matrix(2, 1) = 5;
    matrix(2, 2) = 6;

    REQUIRE_THAT(matrix.toArma(), IsApproxEqual(arma::mat{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}}, 0));
}