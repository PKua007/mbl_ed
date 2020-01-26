//
// Created by Piotr Kubala on 21/01/2020.
//

#ifndef MBL_ED_FRONTEND_H
#define MBL_ED_FRONTEND_H


#include "analyzer/Analyzer.h"
#include "Parameters.h"

class Frontend {
private:
    std::ostream &out;

    auto buildHamiltonianGenerator(const Parameters &params, bool changePhi0ForAverage);
    Analyzer prepareAnalyzer(const std::vector<std::string> &tasks);

    template<template<typename> typename AveragingModel_t, typename HamiltonianGenerator_t>
    void perform_simulations(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, Analyzer &analyzer,
                             std::size_t from, std::size_t to, std::size_t totalSimulations,
                             const std::string &fileSignature, bool saveEigenenergies);

public:
    explicit Frontend(std::ostream &out) : out{out} { }

    void simulate(int argc, char **argv);
    void analyze(int argc, char **argv);
    void printGeneralHelp(const std::string &cmd);
};


#endif //MBL_ED_FRONTEND_H
