//
// Created by Piotr Kubala on 23/09/2020.
//

#include <sstream>

#include "ObservablesBuilder.h"
#include "core/observables/OnsiteOccupations.h"
#include "core/observables/OnsiteOccupationsSquared.h"
#include "core/observables/OnsiteFluctuations.h"
#include "core/observables/Correlations.h"
#include "core/observables/BipariteEntropy.h"
#include "core/observables/CavityOnsiteOccupations.h"
#include "core/observables/CavityOnsiteOccupationsSquared.h"

void ObservablesBuilder::build(const std::vector<std::string> &observables, const Parameters &params,
                               const std::shared_ptr<FockBasis> &fockBasis,
                               std::optional<std::reference_wrapper<const HamiltonianGenerator>> hamiltonianGenerator)
{
    this->primaryObservables.clear();
    this->secondaryObservables.clear();
    this->storedObservables.clear();

    // Here are of all possible PrimaryObservables, but will be initialized only if needed
    std::shared_ptr<OnsiteOccupations> occupations;
    std::shared_ptr<OnsiteOccupationsSquared> occupations2;
    std::shared_ptr<BipariteEntropy> bipariteEntropy;
    std::shared_ptr<CavityOnsiteOccupations> cavityOccupations;
    std::shared_ptr<CavityOnsiteOccupationsSquared> cavityOccupationsSquared;

    for (const auto &observable : observables) {
        std::istringstream observableStream(observable);
        std::string observableName;
        observableStream >> observableName;

        // Here PrimaryObservables are created on explicit request, or if SecondaryObservables need them. They
        // are created only once
        if (observableName == "n_i") {
            occupations = std::make_shared<OnsiteOccupations>(fockBasis);
            this->storedObservables.push_back(occupations);
        } else if (observableName == "n_iN_j") {
            occupations2 = std::make_shared<OnsiteOccupationsSquared>(fockBasis);
            this->storedObservables.push_back(occupations2);
        } else if (observableName == "rho_i") {
            if (occupations == nullptr)
                occupations = std::make_shared<OnsiteOccupations>(fockBasis);
            if (occupations2 == nullptr)
                occupations2 = std::make_shared<OnsiteOccupationsSquared>(fockBasis);
            auto fluctuations = std::make_shared<OnsiteFluctuations>(params.K);
            this->secondaryObservables.push_back(fluctuations);
            this->storedObservables.push_back(fluctuations);
        } else if (observableName == "G_d") {
            std::size_t marginSize{};
            observableStream >> marginSize;
            ValidateMsg(observableStream, "Wrong G_d format. Usage: G_d [margin size]");
            Validate(2*marginSize + 2 <= params.K);
            if (occupations == nullptr)
                occupations = std::make_shared<OnsiteOccupations>(fockBasis);
            if (occupations2 == nullptr)
                occupations2 = std::make_shared<OnsiteOccupationsSquared>(fockBasis);
            auto correlations = std::make_shared<Correlations>(params.K, marginSize);
            this->secondaryObservables.push_back(correlations);
            this->storedObservables.push_back(correlations);
        } else if (observableName == "S") {
            bipariteEntropy = std::make_shared<BipariteEntropy>(fockBasis);
            this->storedObservables.push_back(bipariteEntropy);
        } else if (observableName == "n_i_cos") {
            if (!hamiltonianGenerator.has_value())
                throw std::runtime_error("HamiltonianGenerator is needed to create CavityOccupations");

            std::shared_ptr<CavityLongInteraction> cavityLongInteraction;
            for (const auto &term : hamiltonianGenerator->get().getDiagonalTerms()) {
                cavityLongInteraction = std::dynamic_pointer_cast<CavityLongInteraction>(term);
                if (cavityLongInteraction)
                    break;
            }
            if (!cavityLongInteraction)
                throw std::runtime_error("CavityLongInteraction term not found for CavityOccupations");

            cavityOccupations = std::make_shared<CavityOnsiteOccupations>(fockBasis, cavityLongInteraction);
            this->storedObservables.push_back(cavityOccupations);
        } else if (observableName == "n_iN_j_cos") {
            if (!hamiltonianGenerator.has_value())
                throw std::runtime_error("HamiltonianGenerator is needed to create CavityOccupationsSquared");

            std::shared_ptr<CavityLongInteraction> cavityLongInteraction;
            for (const auto &term : hamiltonianGenerator->get().getDiagonalTerms()) {
                cavityLongInteraction = std::dynamic_pointer_cast<CavityLongInteraction>(term);
                if (cavityLongInteraction)
                    break;
            }
            if (!cavityLongInteraction)
                throw std::runtime_error("CavityLongInteraction term not found for CavityOccupationsSquared");

            cavityOccupationsSquared = std::make_shared<CavityOnsiteOccupationsSquared>(
                fockBasis, cavityLongInteraction
            );
            this->storedObservables.push_back(cavityOccupationsSquared);
        } else {
            throw ValidationException("Unknown observable: " + observable);
        }

        if (occupations != nullptr)
            this->primaryObservables.push_back(occupations);
        if (occupations2 != nullptr)
            this->primaryObservables.push_back(occupations2);
        if (bipariteEntropy != nullptr)
            this->primaryObservables.push_back(bipariteEntropy);
        if (cavityOccupations != nullptr)
            this->primaryObservables.push_back(cavityOccupations);
        if (cavityOccupationsSquared != nullptr)
            this->primaryObservables.push_back(cavityOccupationsSquared);
    }
}
