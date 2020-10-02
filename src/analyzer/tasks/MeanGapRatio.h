//
// Created by pkua on 08.11.2019.
//

#ifndef MBL_ED_MEANGAPRATIO_H
#define MBL_ED_MEANGAPRATIO_H


#include <vector>

#include "core/Eigensystem.h"
#include "utils/Quantity.h"
#include "analyzer/InlineAnalyzerTask.h"

/**
 * @brief InlineAnalyzerTask which computer mean gap ratio of eigenenergies normalized to [0, 1] from a specific energy
 * band.
 *
 * @p Gap ratio \f$ \bar{r} \f$ is the ratio of the gaps between consecutive eigenenergies
 * \f$ \epsilon_1, \dots, \epsilon_\text{dim} \f$. Is is defined to be always non-greater than 1:
 * \f[
 *     \bar{r}(\epsilon_i) = \min
 *         \left(
 *             \frac{\epsilon_{i+1} - \epsilon_i}{\epsilon_i - \epsilon_{i-1}},
 *             \frac{\epsilon_i - \epsilon_{i-1}}{\epsilon_{i+1} - \epsilon_i}
 *         \right)
 * \f]
 *
 * @p Averaging is done in a following way: for a single Eigensystem the average value is computed, and this average
 * is a single sample value in the next averaging performed over multiple Eigensystems. The error is also of the
 * second average.
 */
class MeanGapRatio : public InlineAnalyzerTask {
private:
    std::optional<FockBasis::Vector> middleVector;
    double relativeMiddleEnergy{};
    double relativeMargin{};
    std::vector<double> gapRatios{};

    double calculateEnergyOfFockState(const FockBasis::Vector &state, const Eigensystem &eigensystem) const;
    [[nodiscard]] Quantity calculateMean() const;

public:
    /**
     * @brief Constructs the class, which will compute mean gap ratio only for normalized eigenenergies from a specific
     * energy band.
     * @param relativeMiddleEnergy the middle of the band (in the [0, 1] regime)
     * @param relativeMargin the width of the band (also in the [0, 1] regime)
     */
    MeanGapRatio(double relativeMiddleEnergy, double relativeMargin);

    /**
     * @brief Constructs the class, which will compute mean gap ratio only for normalized eigenerergies around a
     * specific vector.
     * @param middleVector vector around which we should compute mgr
     * @param relativeMargin relativeMargin the width of the band (in the [0, 1] regime)
     */
    MeanGapRatio(const FockBasis::Vector &middleVector, double relativeMargin);

    /**
     * @brief Adds value for a given @a eigensystem to the average mean gap ratio.
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


#endif //MBL_ED_MEANGAPRATIO_H
