//
// Created by pkua on 08.06.2020.
//

#include "CavityConstantsAveragingModel.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/RND.h"
#include "simulation/terms/LookupCavityZ2.h"
#include "simulation/terms/OnsiteDisorder.h"
#include "simulation/terms/LookupCavityYZ.h"
#include "simulation/terms/LookupCavityY2.h"

void CavityConstantsAveragingModel::setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                                              unsigned long simulationIndex,
                                                              unsigned long numberOfSimulations)
{
    Expects(numberOfSimulations > 0);
    Expects(simulationIndex < numberOfSimulations);

    bool termFound{};
    for (auto &diagonalTerm : hamiltonianGenerator.getDiagonalTerms()) {
        try {
            auto &lookupCavityZ2 = dynamic_cast<LookupCavityZ2&>(*diagonalTerm);
            lookupCavityZ2.changeRealisation(simulationIndex);
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }

        try {
            auto &onsiteDisorder = dynamic_cast<OnsiteDisorder&>(*diagonalTerm);
            onsiteDisorder.resampleOnsiteEnergies(rnd);
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }
    }

    for (auto &hoppingTerm : hamiltonianGenerator.getHoppingTerms()) {
        try {
            auto &lookupCavityYZ = dynamic_cast<LookupCavityYZ&>(*hoppingTerm);
            lookupCavityYZ.changeRealisation(simulationIndex);
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }
    }

    for (auto &doubleHoppingTerm : hamiltonianGenerator.getDoubleHoppingTerms()) {
        try {
            auto &lookupCavityY2 = dynamic_cast<LookupCavityY2&>(*doubleHoppingTerm);
            lookupCavityY2.changeRealisation(simulationIndex);
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }
    }

    if (!termFound)
        throw std::runtime_error("CavityConstantsAveragingModel: lack of any term to average on");
}