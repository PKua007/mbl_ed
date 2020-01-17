//
// Created by Piotr Kubala on 16/01/2020.
//

#ifndef MBL_ED_ANALYZERTASK_H
#define MBL_ED_ANALYZERTASK_H

#include <vector>
#include <string>

class AnalyzerTask {
public:
    virtual ~AnalyzerTask() = default;

    virtual void analyze(const std::vector<double> &eigenenergies) = 0;
    [[nodiscard]] virtual std::string getName() const = 0;
};


#endif //MBL_ED_ANALYZERTASK_H