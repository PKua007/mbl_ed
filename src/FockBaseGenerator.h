//
// Created by pkua on 01.11.2019.
//

#ifndef CAVITY_MBL_FOCKBASEGENERATOR_H
#define CAVITY_MBL_FOCKBASEGENERATOR_H

#include <vector>

#include "FockBase.h"

class FockBaseGenerator {
public:
    [[nodiscard]] FockBase generate(int numberOfSites, int numberOfParticles) const;
};


#endif //CAVITY_MBL_FOCKBASEGENERATOR_H
