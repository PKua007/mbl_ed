//
// Created by Piotr Kubala on 22/01/2020.
//

#ifndef MBL_ED_MEANINVERSEPARTICIPATIONRATIO_H
#define MBL_ED_MEANINVERSEPARTICIPATIONRATIO_H

#include "analyzer/InlineAnalyzerTask.h"
#include "analyzer/BandExtractor.h"
#include "utils/Quantity.h"

/**
 * @brief InlineAnalyzerTask which computer inverse participation ratio of eigenenergies normalized to [0, 1] from a
 * specific energy band.
 *
 * @p For the definition of inverse participation ratio see InverseParticipationRatio.
 *
 * @p Averaging is done in a following way: for a single Eigensystem the average value is computed, and this average
 * is a single sample value in the next averaging performed over multiple Eigensystems. The error is also of the
 * second average.
 */
class MeanInverseParticipationRatio : public InlineAnalyzerTask {
private:
    BandExtractor extractor;
    std::vector<double> ratios{};

    [[nodiscard]] Quantity calculateMean() const;

public:
    /**
     * @brief Constructs the class, which will compute inverse participation ratio only for normalized eigenenergies
     * from a specific energy band.
     * @param range the range to choose eigenstates from
     */
    explicit MeanInverseParticipationRatio(BandExtractor::Range range) : extractor(std::move(range), "Mean ipr") { }

    /**
     * @brief Adds a value for a given @a eigensystem to the average inverse participation ratio.
     * @param eigensystem Eigensystem to fetch data points from
     * @param logger unused
     */
    void analyze(const Eigensystem &eigensystem, Logger &logger) override;
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::vector<std::string> getResultHeader() const override;
    [[nodiscard]] std::vector<std::string> getResultFields() const override;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_MEANINVERSEPARTICIPATIONRATIO_H
