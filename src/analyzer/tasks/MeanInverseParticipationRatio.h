//
// Created by Piotr Kubala on 22/01/2020.
//

#ifndef MBL_ED_MEANINVERSEPARTICIPATIONRATIO_H
#define MBL_ED_MEANINVERSEPARTICIPATIONRATIO_H

#include "analyzer/InlineAnalyzerTask.h"
#include "utils/Quantity.h"

/**
 * @brief InlineAnalyzerTask which computer inverse participation ratio of eigenenergies normalized to [0, 1] from a
 * specific energy band.
 *
 * For the definition of inverse participation ratio see InverseParticipationRatio.
 */
class MeanInverseParticipationRatio : public InlineAnalyzerTask {
private:
    double relativeMiddleEnergy{};
    double relativeMargin{};
    std::vector<double> ratios{};

    [[nodiscard]] Quantity calculateMean() const;

public:
    /**
     * @brief Constructs the class, which will compute inverse participation ratio only for normalized eigenenergies
     * from a specific energy band.
     * @param relativeMiddleEnergy the middle of the band (in the [0, 1] regime)
     * @param relativeMargin the width of the band (also in the [0, 1] regime)
     */
    MeanInverseParticipationRatio(double relativeMiddleEnergy, double relativeMargin);

    /**
     * @brief Adds point from a given @a eigensystem to the average inverse participation ratio.
     *
     * Multiple invocations result in adding more points to the average.
     * @param eigensystem Eigensystem to fetch data points from
     * @param logger unused
     */
    void analyze(const Eigensystem &eigensystem, std::ostream &logger) override;
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::vector<std::string> getResultHeader() const override;
    [[nodiscard]] std::vector<std::string> getResultFields() const override;
};


#endif //MBL_ED_MEANINVERSEPARTICIPATIONRATIO_H
