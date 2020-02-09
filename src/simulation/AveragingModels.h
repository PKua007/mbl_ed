//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_AVERAGINGMODELS_H
#define MBL_ED_AVERAGINGMODELS_H

#include <cmath>

#include "simulation/terms/CavityLongInteraction.h"
#include "simulation/terms/OnsiteDisorder.h"
#include "utils/Assertions.h"
#include "HamiltonianGenerator.h"
#include "RND.h"

/**
 * @brief Averaging model, used in Simulation, which samples new onsite disorder terms for each simulation.
 * @tparam DisorderGenerator_t disorder generator used
 */
template<typename DisorderGenerator_t>
class OnsiteDisorderAveragingModel {
public:
    static void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        static_cast<void>(simulationIndex);
        static_cast<void>(numberOfSimulations);

        for (auto &diagonalTerm : hamiltonianGenerator.getDiagonalTerms()) {
            try {
                auto &onsiteDisorder = dynamic_cast<OnsiteDisorder<DisorderGenerator_t>&>(*diagonalTerm);
                onsiteDisorder.resampleOnsiteEnergies(rnd);
            } catch (std::bad_cast &e) { static_cast<void>(e); }
        }
    }
};

/**
 * @brief Averaging model, used in Simulation, which distributes evenly phi0 values across simulations based on
 * @a simulationIndex and @a numberOfSimulations.
 * @details Onsite disorder is also resampled.
 * @tparam DisorderGenerator_t disorder generator used
 */
template<typename DisorderGenerator_t>
class UniformPhi0AveragingModel {
public:
    static void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        Expects(numberOfSimulations > 0);
        Expects(simulationIndex < numberOfSimulations);

        for (auto &diagonalTerm : hamiltonianGenerator.getDiagonalTerms()) {
            try {
                auto &longIteraction = dynamic_cast<CavityLongInteraction&>(*diagonalTerm);
                double phi0 = M_PI * simulationIndex / numberOfSimulations;
                longIteraction.setPhi0(phi0);
            } catch (std::bad_cast &e) { static_cast<void>(e); }

            try {
                auto &onsiteDisorder = dynamic_cast<OnsiteDisorder<DisorderGenerator_t>&>(*diagonalTerm);
                onsiteDisorder.resampleOnsiteEnergies(rnd);
            } catch (std::bad_cast &e) { static_cast<void>(e); }
        }
    }
};

/**
 * @brief Averaging model, used in Simulation, which picks random phi0 values.
 * @details Onsite disorder is also resampled.
 * @tparam DisorderGenerator_t disorder generator used
 */
template<typename DisorderGenerator_t>
class RandomPhi0AveragingModel {
public:
    static void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        Expects(numberOfSimulations > 0);
        Expects(simulationIndex < numberOfSimulations);

        for (auto &diagonalTerm : hamiltonianGenerator.getDiagonalTerms()) {
            try {
                auto &longIteraction = dynamic_cast<CavityLongInteraction&>(*diagonalTerm);
                longIteraction.setPhi0(2*M_PI*rnd());
            } catch (std::bad_cast &e) { static_cast<void>(e); }

            try {
                auto &onsiteDisorder = dynamic_cast<OnsiteDisorder<DisorderGenerator_t>&>(*diagonalTerm);
                onsiteDisorder.resampleOnsiteEnergies(rnd);
            } catch (std::bad_cast &e) { static_cast<void>(e); }
        }
    }
};

#endif //MBL_ED_AVERAGINGMODELS_H
