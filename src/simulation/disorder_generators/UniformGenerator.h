//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_UNIFORMGENERATOR_H
#define MBL_ED_UNIFORMGENERATOR_H


#include "simulation/DisorderGenerator.h"

/**
 * @brief Disorder generator wich just samples random number from [min, max).
 */
class UniformGenerator : public DisorderGenerator {
private:
    double min{};
    double max{};

public:
    UniformGenerator(double min, double max) : min{min}, max{max} { }

    double generate(RND &rnd) override {
        return this->min + rnd()*(this->max - this->min);
    }
};



#endif //MBL_ED_UNIFORMGENERATOR_H
