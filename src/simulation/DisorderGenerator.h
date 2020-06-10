//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_DISORDERGENERATOR_H
#define MBL_ED_DISORDERGENERATOR_H

#include "RND.h"

/**
 * @brief A class sampling random chemical potential for OnsiteDisorder class from a giver distribution.
 */
class DisorderGenerator {
public:
    /**
     * @brief Sample next chemical potential using passed @a rnd.
     */
    virtual double generate(RND &rnd) = 0;
};

#endif //MBL_ED_DISORDERGENERATOR_H
