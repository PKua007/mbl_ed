//
// Created by pkua on 08.06.2020.
//

#ifndef MBL_ED_ANALYZERBUILDER_H
#define MBL_ED_ANALYZERBUILDER_H

#include <vector>

#include "analyzer/Analyzer.h"
#include "Parameters.h"

class AnalyzerBuilder {
public:
    Analyzer build(const std::vector<std::string> &tasks, const Parameters &params,
                                    std::shared_ptr<FockBase> fockBase);
};


#endif //MBL_ED_ANALYZERBUILDER_H
