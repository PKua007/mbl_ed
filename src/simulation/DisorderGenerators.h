//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_DISORDERGENERATORS_H
#define MBL_ED_DISORDERGENERATORS_H

#include "RND.h"

/**
 * @brief Disorder generator wich just samples random number from [min, max).
 */
class UniformGenerator {
private:
    double min{};
    double max{};

public:
    UniformGenerator(double min, double max) : min{min}, max{max} { }

    double operator()(RND &rnd) {
        return this->min + rnd()*(this->max - this->min);
    }
};

#endif //MBL_ED_DISORDERGENERATORS_H
