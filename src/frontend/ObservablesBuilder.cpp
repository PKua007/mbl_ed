//
// Created by Piotr Kubala on 23/09/2020.
//

#include <sstream>

#include "ObservablesBuilder.h"
#include "evolution/observables/OnsiteOccupations.h"
#include "evolution/observables/OnsiteOccupationsSquared.h"
#include "evolution/observables/OnsiteFluctuations.h"
#include "evolution/observables/Correlations.h"

void ObservablesBuilder::build(const std::vector<std::string> &observables, const Parameters &params,
                               const std::shared_ptr<FockBase> &fockBase)
{
    this->primaryObservables.clear();
    this->secondaryObservables.clear();
    this->storedObservables.clear();

    std::shared_ptr<OnsiteOccupations> occupations;
    std::shared_ptr<OnsiteOccupationsSquared> occupations2;

    for (const auto &observable : observables) {
        std::istringstream observableStream(observable);
        std::string observableName;
        observableStream >> observableName;

        if (observableName == "n_i") {
            occupations = std::make_shared<OnsiteOccupations>(params.K, fockBase);
            this->storedObservables.push_back(occupations);
        } else if (observableName == "n_iN_j") {
            occupations2 = std::make_shared<OnsiteOccupationsSquared>(params.K, fockBase);
            this->storedObservables.push_back(occupations2);
        } else if (observableName == "rho_i") {
            if (occupations == nullptr)
                occupations = std::make_shared<OnsiteOccupations>(params.K, fockBase);
            if (occupations2 == nullptr)
                occupations2 = std::make_shared<OnsiteOccupationsSquared>(params.K, fockBase);
            auto fluctuations = std::make_shared<OnsiteFluctuations>(params.K);
            this->secondaryObservables.push_back(fluctuations);
            this->storedObservables.push_back(fluctuations);
        } else if (observableName == "G_d") {
            std::size_t marginSize{};
            observableStream >> marginSize;
            ValidateMsg(observableStream, "Wrong G_d format. Usage: G_d [margin size]");
            Validate(2*marginSize + 2 <= params.K);
            if (occupations == nullptr)
                occupations = std::make_shared<OnsiteOccupations>(params.K, fockBase);
            if (occupations2 == nullptr)
                occupations2 = std::make_shared<OnsiteOccupationsSquared>(params.K, fockBase);
            auto correlations = std::make_shared<Correlations>(params.K, marginSize);
            this->secondaryObservables.push_back(correlations);
            this->storedObservables.push_back(correlations);
        } else {
            throw ValidationException("Unknown observable: " + observable);
        }

        if (occupations != nullptr)
            this->primaryObservables.push_back(occupations);
        if (occupations2 != nullptr)
            this->primaryObservables.push_back(occupations2);
    }
}
