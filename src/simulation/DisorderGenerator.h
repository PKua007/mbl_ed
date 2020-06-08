//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_DISORDERGENERATOR_H
#define MBL_ED_DISORDERGENERATOR_H

#include "RND.h"

class DisorderGenerator {
public:
    virtual double generate(RND &rnd) = 0;
};

#endif //MBL_ED_DISORDERGENERATOR_H
