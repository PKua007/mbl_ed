//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_INLINEANALYZERTASK_H
#define MBL_ED_INLINEANALYZERTASK_H


#include "AnalyzerTask.h"

class InlineAnalyzerTask : public AnalyzerTask {
public:
    [[nodiscard]] virtual std::vector<std::string> getResultHeader() const = 0;
    [[nodiscard]] virtual std::vector<std::string> getResultFields() const = 0;
};


#endif //MBL_ED_INLINEANALYZERTASK_H
