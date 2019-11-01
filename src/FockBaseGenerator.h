//
// Created by pkua on 01.11.2019.
//

#ifndef CAVITY_MBL_FOCKBASEGENERATOR_H
#define CAVITY_MBL_FOCKBASEGENERATOR_H

#include <vector>

class FockBaseGenerator {
public:
    std::vector<std::vector<int>> generate(int M, int N) const;
};


#endif //CAVITY_MBL_FOCKBASEGENERATOR_H
