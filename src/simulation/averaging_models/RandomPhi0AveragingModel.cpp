//
// Created by pkua on 08.06.2020.
//

#include "RandomPhi0AveragingModel.h"

void RandomPhi0AveragingModel::setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                                         std::size_t simulationIndex, std::size_t numberOfSimulations)
{
    Expects(numberOfSimulations > 0);
    Expects(simulationIndex < numberOfSimulations);

    bool termFound{};
    for (auto &diagonalTerm : hamiltonianGenerator.getDiagonalTerms()) {
        try {
            auto &quasiperiodicDisorder = dynamic_cast<QuasiperiodicDisorder&>(*diagonalTerm);
            quasiperiodicDisorder.setPhi0(2*M_PI*rnd());
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }
        try {
            auto &longIteraction = dynamic_cast<CavityLongInteraction&>(*diagonalTerm);
            longIteraction.setPhi0(2*M_PI*rnd());
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }

        try {
            auto &onsiteDisorder = dynamic_cast<OnsiteDisorder&>(*diagonalTerm);
            onsiteDisorder.resampleOnsiteEnergies(rnd);
            termFound = true;
        } catch (std::bad_cast &e) { static_cast<void>(e); }
    }

    if (!termFound)
        throw std::runtime_error("RandomPhi0AveragingModel: lack of any term to average on");
}
