//
// Created by pkua on 08.11.2019.
//

#ifndef MBL_ED_MEANGAPRATIO_H
#define MBL_ED_MEANGAPRATIO_H


#include <utility>
#include <vector>
#include <variant>

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
public:
    struct EpsilonRange {
        EpsilonRange() = default;
        EpsilonRange(double epsilonMiddle, double epsilonMargin);

        double epsilonMiddle{};
        double epsilonMargin{};
    };

    struct VectorRange {
        VectorRange() = default;
        VectorRange(FockBasis::Vector middleVector, double epsilonMargin);

        FockBasis::Vector middleVector{};
        double epsilonMargin{};
    };

    struct CDFRange {
        CDFRange() = default;
        CDFRange(double cdfMiddle, double cdfMargin);

        double cdfMiddle{};
        double cdfMargin{};
    };

    using Range = std::variant<EpsilonRange, VectorRange, CDFRange>;

private:
    Range range;
    std::vector<double> gapRatios{};

    [[nodiscard]] double calculateEnergyOfFockState(const FockBasis::Vector &state,
                                                    const Eigensystem &eigensystem) const;
    [[nodiscard]] Quantity calculateMean() const;

public:
    explicit MeanGapRatio(Range range) : range{std::move(range)} { }

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
