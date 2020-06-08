//
// Created by Piotr Kubala on 21/01/2020.
//

#ifndef MBL_ED_FRONTEND_H
#define MBL_ED_FRONTEND_H


#include <evolution/CorrelationsTimeEvolutionParameters.h>
#include "analyzer/Analyzer.h"
#include "Parameters.h"
#include "simulation/SimulationParameters.h"
#include "simulation/RND.h"
#include "simulation/HamiltonianGenerator.h"

/**
 * @brief Class responsible for communication between the user and simulation backend.
 */
class Frontend {
private:
    std::ostream &out;

    Analyzer prepareAnalyzer(const std::vector<std::string> &tasks, const Parameters &params,
                             std::shared_ptr<FockBase> fockBase);

    template<template <typename> typename AveragingModel_t>
    void perform_simulations(std::unique_ptr<HamiltonianGenerator> hamiltonianGenerator, std::unique_ptr<RND>,
                             Analyzer &analyzer, const SimulationParameters &simulationParameters);

    template<template <typename> typename AveragingModel_t>
    void perform_chebyshev_evolution(std::unique_ptr<HamiltonianGenerator> hamiltonianGenerator,
                                     std::unique_ptr<RND> rnd, const Parameters &params,
                                     const CorrelationsTimeEvolutionParameters &evolutionParameters);

public:
    explicit Frontend(std::ostream &out) : out{out} { }

    void simulate(int argc, char **argv);
    void analyze(int argc, char **argv);
    void chebyshev(int argc, char **argv);
    void printGeneralHelp(const std::string &cmd);
};


#endif //MBL_ED_FRONTEND_H
