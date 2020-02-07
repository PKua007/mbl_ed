//
// Created by Piotr Kubala on 07/02/2020.
//

#ifndef MBL_ED_AVERAGINGMODELS_H
#define MBL_ED_AVERAGINGMODELS_H

#include <cmath>

#include "utils/Assertions.h"

template<typename HamiltonianGenerator_t>
class OnsiteDisorderAveragingModel {
public:
    static void setupHamiltonianGenerator(HamiltonianGenerator_t &hamiltonianGenerator, std::size_t simulationIndex,
                                          std::size_t numberOfSimulations)
    {
        static_cast<void>(simulationIndex);
        static_cast<void>(numberOfSimulations);
        hamiltonianGenerator.resampleOnsiteEnergies();
    }
};

template<typename HamiltonianGenerator_t>
class Phi0AveragingModel {
public:
    static void setupHamiltonianGenerator(HamiltonianGenerator_t &hamiltonianGenerator,
                                          std::size_t simulationIndex, std::size_t numberOfSimulations)
    {
        Expects(numberOfSimulations > 0);
        Expects(simulationIndex < numberOfSimulations);

        double phi0 = 2*M_PI*simulationIndex/numberOfSimulations;
        hamiltonianGenerator.resampleOnsiteEnergies();
        hamiltonianGenerator.setPhi0(phi0);
    }
};

#endif //MBL_ED_AVERAGINGMODELS_H
