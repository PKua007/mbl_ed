//
// Created by pkua on 08.11.2019.
//

#ifndef MBL_ED_MEANGAPRATIO_H
#define MBL_ED_MEANGAPRATIO_H


#include <vector>
#include <simulation/Eigensystem.h>

#include "utils/Quantity.h"
#include "analyzer/InlineAnalyzerTask.h"

class MeanGapRatio : public InlineAnalyzerTask {
private:
    double relativeMiddleEnergy{};
    double relativeMargin{};
    std::vector<double> gapRatios{};

    [[nodiscard]] std::vector<double> getNormalizedEigenenergies(const std::vector<double> &eigenenergies) const;
    [[nodiscard]] Quantity calculateMean() const;

public:
    MeanGapRatio(double relativeMiddleEnergy, double relativeMargin);

    void analyze(const Eigensystem &eigensystem) override;
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::vector<std::string> getResultHeader() const override;
    [[nodiscard]] std::vector<std::string> getResultFields() const override;
};


#endif //MBL_ED_MEANGAPRATIO_H
