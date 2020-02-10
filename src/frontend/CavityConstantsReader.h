//
// Created by Piotr Kubala on 10/02/2020.
//

#ifndef MBL_ED_CAVITYCONSTANTSREADER_H
#define MBL_ED_CAVITYCONSTANTSREADER_H

#include <istream>

#include "simulation/CavityConstants.h"

class CavityConstantsReader {
public:
    static CavityConstants load(std::istream &in);
};


#endif //MBL_ED_CAVITYCONSTANTSREADER_H
