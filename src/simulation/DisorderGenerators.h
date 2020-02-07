//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_DISORDERGENERATORS_H
#define MBL_ED_DISORDERGENERATORS_H

#include <random>

/**
 * @brief Disorder generator wich just samples random number from [min, max).
 */
class UniformGenerator {
private:
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;

public:
    UniformGenerator(double min, double max, unsigned long seed) : generator(seed), distribution(min, max) { }

    UniformGenerator(UniformGenerator &dummy) = delete;
    UniformGenerator operator=(UniformGenerator dummy) = delete;

    double operator()() {
        return this->distribution(this->generator);
    }
};

#endif //MBL_ED_DISORDERGENERATORS_H
