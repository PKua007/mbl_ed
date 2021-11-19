//
// Created by Piotr Kubala on 21/01/2020.
//

#ifndef MBL_ED_FRONTEND_H
#define MBL_ED_FRONTEND_H


#include "evolution/TimeEvolutionParameters.h"
#include "analyzer/Analyzer.h"
#include "Parameters.h"
#include "simulation/ExactDiagonalizationParameters.h"
#include "core/RND.h"
#include "core/HamiltonianGenerator.h"

/**
 * @brief Class responsible for the communication between the user and the simulation backend.
 */
class Frontend {
private:
    std::ostream &out;

    void setOverridenParamsAsAdditionalText(Logger &logger, std::vector<std::string> overridenParams) const;
    void setVerbosityLevel(Logger &logger, const std::string &verbosityLevelName) const;

    [[nodiscard]] ExactDiagonalizationParameters
    prepareExactDiagonalizationParameters(const std::filesystem::path &directory, const Parameters &params) const;

    [[nodiscard]] std::string stripSuffix(const std::string &string, const std::string &suffix) const;

    void appendNotMatchingSignatures(std::vector<std::string> &notMatchingSignatures,
                                     const std::vector<std::string> &filenames1, const std::string &suffix1,
                                     const std::vector<std::string> &filenames2, const std::string &suffix2);

    void validateEigensystemFiles(Logger &logger, const std::vector<std::string> &energiesFilenames,
                                  const std::vector<std::string> &statesFilenames);

public:
    explicit Frontend(std::ostream &out) : out{out} { }

    void ed(int argc, char **argv);
    void analyze(int argc, char **argv);
    void chebyshev(int argc, char **argv);
    void quench(int argc, char **argv);
    void randomStates(int argc, char **argv);
    void preview(int argc, char **argv);

    void printGeneralHelp(const std::string &cmd);

    Eigensystem
    restoreEigensystem(const std::string &energiesFilename, bool restoreEigenstates, std::shared_ptr<FockBasis> basis,
                       Logger &logger) const;
};


#endif //MBL_ED_FRONTEND_H
