//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_ANALYZERBUILDER_H
#define MBL_ED_ANALYZERBUILDER_H

#include <vector>

#include "analyzer/Analyzer.h"
#include "Parameters.h"

/**
 * @brief A class responsible for building the Analyzer class based on the vector of given tasks and Parameters.
 * @details Extracted from Frontend mainly to reduce compilation times.
 */
class AnalyzerBuilder {
public:
    /**
     * @brief Prepares the Analyzer. Parsing @a tasks includes also validation and throwing exceptions if needed.
     */
    std::unique_ptr<Analyzer> build(const std::vector<std::string> &tasks, const Parameters &params,
                                    std::shared_ptr<FockBasis> fockBasis);
};


#endif //MBL_ED_ANALYZERBUILDER_H
