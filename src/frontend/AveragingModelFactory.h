//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_AVERAGINGMODELFACTORY_H
#define MBL_ED_AVERAGINGMODELFACTORY_H

#include <memory>
#include <string>

#include "simulation/AveragingModel.h"

class AveragingModelFactory {
public:
    std::unique_ptr<AveragingModel> create(const std::string &name);
};


#endif //MBL_ED_AVERAGINGMODELFACTORY_H
