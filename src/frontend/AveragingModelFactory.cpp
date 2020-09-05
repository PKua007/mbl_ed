//
// Created by pkua on 08.06.2020.
//

#include "AveragingModelFactory.h"

#include "core/averaging_models/DummyAveragingModel.h"
#include "core/averaging_models/UniformPhi0AveragingModel.h"
#include "core/averaging_models/RandomPhi0AveragingModel.h"
#include "core/averaging_models/OnsiteDisorderAveragingModel.h"
#include "core/averaging_models/CavityConstantsAveragingModel.h"

std::unique_ptr<AveragingModel> AveragingModelFactory::create(const std::string &name) {
    if (name == "none")
        return std::make_unique<DummyAveragingModel>();
    else if (name == "uniformPhi0")
        return std::make_unique<UniformPhi0AveragingModel>();
    else if (name == "randomPhi0")
        return std::make_unique<RandomPhi0AveragingModel>();
    else if (name == "onsiteDisorder")
        return std::make_unique<OnsiteDisorderAveragingModel>();
    else if (name == "cavityConstants")
        return std::make_unique<CavityConstantsAveragingModel>();
    else
        throw std::runtime_error("Unknown averaging model: " + name);
}
