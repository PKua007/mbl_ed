//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_AVERAGINGMODELFACTORY_H
#define MBL_ED_AVERAGINGMODELFACTORY_H

#include <memory>
#include <string>

#include "simulation/AveragingModel.h"

/**
 * @brief A class creating appropriate AveragingModel based on its name.
 * @details Extracted from Frontend mainly to decrease compilation times.
 */
class AveragingModelFactory {
public:
    std::unique_ptr<AveragingModel> create(const std::string &name);
};


#endif //MBL_ED_AVERAGINGMODELFACTORY_H
