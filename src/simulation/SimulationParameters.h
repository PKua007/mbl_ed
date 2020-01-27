//
// Created by Piotr Kubala on 27/01/2020.
//

#ifndef MBL_ED_SIMULATIONPARAMETERS_H
#define MBL_ED_SIMULATIONPARAMETERS_H

struct SimulationParameters {
    std::size_t from{};
    std::size_t to{};
    std::size_t totalSimulations{};
    bool calculateEigenvectors{};
    bool saveEigenenergies{};
    std::string fileSignature{};
};

#endif //MBL_ED_SIMULATIONPARAMETERS_H
