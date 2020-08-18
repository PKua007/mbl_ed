//
// Created by Piotr Kubala on 21/01/2020.
//

#ifndef MBL_ED_FRONTEND_H
#define MBL_ED_FRONTEND_H


#include "evolution/CorrelationsTimeEvolutionParameters.h"
#include "analyzer/Analyzer.h"
#include "Parameters.h"
#include "simulation/SimulationParameters.h"
#include "simulation/RND.h"
#include "simulation/HamiltonianGenerator.h"

/**
 * @brief Class responsible for the communication between the user and the simulation backend.
 */
class Frontend {
private:
    std::ostream &out;

public:
    explicit Frontend(std::ostream &out) : out{out} { }

    void simulate(int argc, char **argv);
    void analyze(int argc, char **argv);
    void chebyshev(int argc, char **argv);
    void quench(int argc, char **argv);
    void printGeneralHelp(const std::string &cmd);
};


#endif //MBL_ED_FRONTEND_H
