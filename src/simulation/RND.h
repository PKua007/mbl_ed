//
// Created by Piotr Kubala on 08/02/2020.
//

#ifndef MBL_ED_RND_H
#define MBL_ED_RND_H

#include <random>

/**
 * @brief A simple wrapper of Mersene Twister @a std::mt19937 generator.
 */
class RND {
private:
    std::uniform_real_distribution<double> distribution;
    std::mt19937 randomGenerator;

public:
    RND() : randomGenerator(std::random_device{}()) { }
    explicit RND(unsigned long seed) : randomGenerator(seed) { }
    virtual ~RND() = default;

    /**
     * @brief Returns random number from [0, 1) interval. Can be mocked.
     */
    virtual double getDouble() { return this->distribution(this->randomGenerator); }

    /**
     * @brief Returns random number from [0, 1) interval. Can be mocked by overriding getDouble().
     */
    double operator()() { return this->getDouble(); }
};


#endif //MBL_ED_RND_H
