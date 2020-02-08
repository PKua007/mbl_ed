//
// Created by Piotr Kubala on 08/02/2020.
//

#ifndef MBL_ED_RND_H
#define MBL_ED_RND_H

#include <random>

class RND {
private:
    std::mt19937 randomGenerator;

public:
    RND() : randomGenerator(std::random_device{}()) { }
    explicit RND(unsigned long seed) : randomGenerator(seed) { }

    double operator()() { return this->randomGenerator(); }
};


#endif //MBL_ED_RND_H
