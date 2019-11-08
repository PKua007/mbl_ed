//
// Created by pkua on 08.11.2019.
//

#ifndef MBL_ED_GAPRATIOCALCULATOR_H
#define MBL_ED_GAPRATIOCALCULATOR_H


#include <vector>
#include "Quantity.h"

class GapRatioCalculator {
private:
    double relativeMiddleEnergy{};
    double relativeMargin{};
    std::vector<double> gapRatios{};

public:
    GapRatioCalculator(double relativeMiddleEnergy, double relativeMargin);

    void addEigenenergies(const std::vector<double> &eigenenergies);
    Quantity calculateMean();

    std::vector<double> getNormalizedEigenenergies(const std::vector<double> &eigenenergies) const;
};


#endif //MBL_ED_GAPRATIOCALCULATOR_H
