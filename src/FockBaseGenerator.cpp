//
// Created by pkua on 01.11.2019.
//

#include <numeric>

#include "FockBaseGenerator.h"

FockBase FockBaseGenerator::generate(int M, int N) const {
    FockBase base;

    // An algorithm from https://arxiv.org/pdf/1102.4006.pdf
    std::vector<int> current(M, 0);
    current[0] = N;
    base.add(current);

    while (current.back() != N) {
        int lastNonzeroK = M - 2;
        while (current[lastNonzeroK] == 0)
            lastNonzeroK--;

        current[lastNonzeroK]--;
        current[lastNonzeroK + 1] = N - std::accumulate(current.begin(), current.begin() + lastNonzeroK + 1, 0);
        std::fill(current.begin() + lastNonzeroK + 2, current.end(), 0);

        base.add(current);
    }

    return base;
}
