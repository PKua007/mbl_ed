//
// Created by pkua on 01.11.2019.
//

#include <numeric>

#include "FockBaseGenerator.h"

std::unique_ptr<FockBase> FockBaseGenerator::generate(int numberOfSites, int numberOfParticles) const {
    auto base = std::make_unique<FockBase>();

    // An algorithm from https://arxiv.org/pdf/1102.4006.pdf
    std::vector<int> current(numberOfSites, 0);
    current[0] = numberOfParticles;
    base->add(current);

    while (current.back() != numberOfParticles) {
        int lastNonzeroK = numberOfSites - 2;
        while (current[lastNonzeroK] == 0)
            lastNonzeroK--;

        current[lastNonzeroK]--;
        current[lastNonzeroK + 1] = numberOfParticles - std::accumulate(current.begin(),
                                                                        current.begin() + lastNonzeroK + 1, 0);
        std::fill(current.begin() + lastNonzeroK + 2, current.end(), 0);

        base->add(current);
    }

    return base;
}
