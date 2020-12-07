//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_ANALYZERBUILDER_H
#define MBL_ED_ANALYZERBUILDER_H

#include <vector>
#include <filesystem>

#include "analyzer/Analyzer.h"
#include "Parameters.h"
#include "core/HamiltonianGenerator.h"

/**
 * @brief A class responsible for building the Analyzer class based on the vector of given tasks and Parameters.
 * @details Extracted from Frontend mainly to reduce compilation times.
 */
class AnalyzerBuilder {
public:
    /**
     * @brief Prepares the Analyzer. Parsing @a tasks includes also validation and throwing exceptions if needed.
     * @detail Optional HamiltonianGenerator may come handy for some tasks, for example using observables relying
     * on specific Hamiltonian terms.
     */
    std::unique_ptr<Analyzer> build(const std::vector<std::string> &tasks, const Parameters &params,
                                    const std::shared_ptr<FockBasis>& fockBasis,
                                    std::optional<std::reference_wrapper<const HamiltonianGenerator>>
                                    hamiltonianGenerator,
                                    const std::filesystem::path &auxiliaryDir = std::filesystem::current_path());
};


#endif //MBL_ED_ANALYZERBUILDER_H
