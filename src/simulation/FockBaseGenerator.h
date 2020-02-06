//
// Created by pkua on 01.11.2019.
//

#ifndef CAVITY_MBL_FOCKBASEGENERATOR_H
#define CAVITY_MBL_FOCKBASEGENERATOR_H

#include <vector>
#include <memory>

#include "FockBase.h"

/**
 * @brief A class generating complete bosonic FockBase of a given number of particles on a given number of sites.
 */
class FockBaseGenerator {
public:
    [[nodiscard]] std::unique_ptr<FockBase> generate(int numberOfSites, int numberOfParticles) const;
};


#endif //CAVITY_MBL_FOCKBASEGENERATOR_H
