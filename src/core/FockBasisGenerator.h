//
// Created by pkua on 01.11.2019.
//

#ifndef FOCK_BASIS_GENERATOR_H
#define FOCK_BASIS_GENERATOR_H

#include <vector>
#include <memory>

#include "FockBasis.h"

/**
 * @brief A class generating complete bosonic FockBasis of a given number of particles on a given number of sites.
 */
class FockBasisGenerator {
public:
    [[nodiscard]] std::unique_ptr<FockBasis> generate(int numberOfParticles, int numberOfSites) const;
};


#endif //FOCK_BASIS_GENERATOR_H
