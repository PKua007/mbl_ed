//
// Created by Piotr Kubala on 22/01/2020.
//

#ifndef MBL_ED_INVERSEPARTICIPATIONRATIO_H
#define MBL_ED_INVERSEPARTICIPATIONRATIO_H

#include "analyzer/InlineAnalyzerTask.h"
#include "utils/Quantity.h"

class InverseParticipationRatio : public InlineAnalyzerTask {
private:
    double relativeMiddleEnergy{};
    double relativeMargin{};
    std::vector<double> ratios{};

    [[nodiscard]] Quantity calculateMean() const;

public:
    InverseParticipationRatio(double relativeMiddleEnergy, double relativeMargin);

    void analyze(const Eigensystem &eigensystem) override;
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::vector<std::string> getResultHeader() const override;
    [[nodiscard]] std::vector<std::string> getResultFields() const override;
};


#endif //MBL_ED_INVERSEPARTICIPATIONRATIO_H
