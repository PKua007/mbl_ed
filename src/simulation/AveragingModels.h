//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_AVERAGINGMODELS_H
#define MBL_ED_AVERAGINGMODELS_H

#include <cmath>
#include <simulation/terms/QuasiperiodicDisorder.h>

#include "simulation/terms/LookupCavityZ2.h"
#include "simulation/terms/LookupCavityYZ.h"
#include "simulation/terms/LookupCavityY2.h"
#include "simulation/terms/CavityLongInteraction.h"
#include "simulation/terms/QuasiperiodicDisorder.h"
#include "simulation/terms/OnsiteDisorder.h"

#include "utils/Assertions.h"
#include "HamiltonianGenerator.h"
#include "RND.h"

/**
 * @brief Averaging model, which does nothing.
 */
template<typename DisorderGenerator_t>
class DummyAveragingModel {
public:
    static void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        static_cast<void>(simulationIndex);
        static_cast<void>(numberOfSimulations);
        static_cast<void>(rnd);
        static_cast<void>(hamiltonianGenerator);
    }
};

/**
 * @brief Averaging model, used in Simulation, which samples new onsite disorder terms for each simulation.
 * @tparam DisorderGenerator_t disorder generator used
 */
template<typename DisorderGenerator_t>
class OnsiteDisorderAveragingModel {
public:
    /**
     * @brief The method looks for OnsiteDisorder term in the @a hamiltonianGenerator.
     */
    static void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        static_cast<void>(simulationIndex);
        static_cast<void>(numberOfSimulations);

        bool changed{};
        for (auto &diagonalTerm : hamiltonianGenerator.getDiagonalTerms()) {
            try {
                auto &onsiteDisorder = dynamic_cast<OnsiteDisorder<DisorderGenerator_t>&>(*diagonalTerm);
                onsiteDisorder.resampleOnsiteEnergies(rnd);
                changed = true;
            } catch (std::bad_cast &e) { static_cast<void>(e); }
        }

        if (!changed)
            throw std::runtime_error("OnsiteDisorderAveragingModel: lack of any term to average on");
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
    /**
     * @brief The method looks for CavityLongInteraction and OnsiteDisorder terms in the @a hamiltonianGenerator.
     */
    static void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        Expects(numberOfSimulations > 0);
        Expects(simulationIndex < numberOfSimulations);

        bool termFound{};
        for (auto &diagonalTerm : hamiltonianGenerator.getDiagonalTerms()) {
            try {
                auto &longIteraction = dynamic_cast<CavityLongInteraction&>(*diagonalTerm);
                double phi0 = M_PI * simulationIndex / numberOfSimulations;
                longIteraction.setPhi0(phi0);
                termFound = true;
            } catch (std::bad_cast &e) { static_cast<void>(e); }

            try {
                auto &quasiperiodicDisorder = dynamic_cast<QuasiperiodicDisorder&>(*diagonalTerm);
                double phi0 = 2 * M_PI * simulationIndex / numberOfSimulations;
                quasiperiodicDisorder.setPhi0(phi0);
                termFound = true;
            } catch (std::bad_cast &e) { static_cast<void>(e); }

            try {
                auto &onsiteDisorder = dynamic_cast<OnsiteDisorder<DisorderGenerator_t>&>(*diagonalTerm);
                onsiteDisorder.resampleOnsiteEnergies(rnd);
                termFound = true;
            } catch (std::bad_cast &e) { static_cast<void>(e); }
        }

        if (!termFound)
            throw std::runtime_error("UniformPhi0AveragingModel: lack of any term to average on");
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
    /**
     * @brief The method looks for CavityLongInteraction and OnsiteDisorder terms in the @a hamiltonianGenerator.
     */
    static void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
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
                auto &onsiteDisorder = dynamic_cast<OnsiteDisorder<DisorderGenerator_t>&>(*diagonalTerm);
                onsiteDisorder.resampleOnsiteEnergies(rnd);
                termFound = true;
            } catch (std::bad_cast &e) { static_cast<void>(e); }
        }

        if (!termFound)
            throw std::runtime_error("RandomPhi0AveragingModel: lack of any term to average on");
    }
};

/**
 * @brief Averaging model, used in Simulation which picks subsequent realisations from CavityConstants.
 * @details It also samples random onsite disorder.
 * @tparam DisorderGenerator_t disorder generator used
 */
template<typename DisorderGenerator_t>
class CavityConstantsAveragingModel {
public:
    /**
     * @brief The method looks for LookupCavityZ2, LookupCavityYZ and OnsiteDisorder terms in the
     * @a hamiltonianGenerator.
     */
    static void setupHamiltonianGenerator(HamiltonianGenerator &hamiltonianGenerator, RND &rnd,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
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
                auto &onsiteDisorder = dynamic_cast<OnsiteDisorder<DisorderGenerator_t>&>(*diagonalTerm);
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
};

#endif //MBL_ED_AVERAGINGMODELS_H
