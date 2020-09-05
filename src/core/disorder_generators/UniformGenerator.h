//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_UNIFORMGENERATOR_H
#define MBL_ED_UNIFORMGENERATOR_H


#include "core/DisorderGenerator.h"

/**
 * @brief DisorderGenerator which just samples uniform number from [@a min, @a max).
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
