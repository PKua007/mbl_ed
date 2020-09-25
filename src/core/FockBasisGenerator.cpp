//
// Created by pkua on 01.11.2019.
//

#include <numeric>

#include "FockBasisGenerator.h"
#include "utils/Assertions.h"

std::unique_ptr<FockBasis> FockBasisGenerator::generate(int numberOfParticles, int numberOfSites) const {
    Expects(numberOfSites > 0);
    auto basis = std::make_unique<FockBasis>();

    // An algorithm from https://arxiv.org/pdf/1102.4006.pdf
    FockBasis::Vector current(numberOfSites, 0);
    current[0] = numberOfParticles;
    basis->add(current);

    while (current.back() != numberOfParticles) {
        int lastNonzeroK = numberOfSites - 2;
        while (current[lastNonzeroK] == 0)
            lastNonzeroK--;

        current[lastNonzeroK]--;
        current[lastNonzeroK + 1] = numberOfParticles - std::accumulate(current.begin(),
                                                                        current.begin() + lastNonzeroK + 1, 0);
        std::fill(current.begin() + lastNonzeroK + 2, current.end(), 0);

        basis->add(current);
    }

    return basis;
}
