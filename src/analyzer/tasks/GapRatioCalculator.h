//
// Created by pkua on 08.11.2019.
//

#ifndef MBL_ED_GAPRATIOCALCULATOR_H
#define MBL_ED_GAPRATIOCALCULATOR_H


#include <vector>

#include "utils/Quantity.h"
#include "analyzer/InlineAnalyzerTask.h"

class GapRatioCalculator : public InlineAnalyzerTask {
private:
    double relativeMiddleEnergy{};
    double relativeMargin{};
    std::vector<double> gapRatios{};

    [[nodiscard]] std::vector<double> getNormalizedEigenenergies(const std::vector<double> &eigenenergies) const;
    [[nodiscard]] Quantity calculateMean() const;

public:
    GapRatioCalculator(double relativeMiddleEnergy, double relativeMargin);

    void analyze(const std::vector<double> &eigenenergies) override;
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::vector<std::string> getResultHeader() const override;
    [[nodiscard]] std::vector<std::string> getResultFields() const override;
};


#endif //MBL_ED_GAPRATIOCALCULATOR_H
