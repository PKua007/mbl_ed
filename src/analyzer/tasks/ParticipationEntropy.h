//
// Created by Piotr Kubala on 06/03/2021.
//

#ifndef MBL_ED_PARTICIPATIONENTROPY_H
#define MBL_ED_PARTICIPATIONENTROPY_H

#include "analyzer/InlineAnalyzerTask.h"
#include "analyzer/BandExtractor.h"
#include "utils/Quantity.h"

/**
 * @brief InlineAnalyzerTask calculating participation entropy.
 *
 * @details <p> For a given eigenvector \f$v^i\f$, \f$S_q\f$ is defined as
 * \f$ \frac{1}{1 - q} \ln \sum_{i=0}^\text{dim}|v^i|^{2q} \f$.
 *
 * <p> For q = 1 it reduces to Shannon's entropy, for q = 2 it reduces to logarithm of InverseParticipationRatio.
 */
class ParticipationEntropy : public InlineAnalyzerTask {
private:
    BandExtractor extractor;
    double q{};
    std::vector<double> entropies{};

    [[nodiscard]] Quantity calculateMean() const;

public:
    /**
     * @brief Constructs the class, which will compute participation entropy only for normalized eigenenergies from a
     * specific energy band.
     * @param relativeMiddleEnergy the middle of the band (in the [0, 1] regime)
     * @param relativeMargin the width of the band (also in the [0, 1] regime)
     */
    ParticipationEntropy(double q, BandExtractor::Range range);

    /**
     * @brief Adds another S_q point to the average from a given @a eigensystem.
     * @details As other similar observables, fist average is over a given energy band from the constructor. Then this
     * is the point for final averaging and calculating final error from many eigensystems.
     */
    void analyze(const Eigensystem &eigensystem, Logger &logger) override;
    [[nodiscard]] std::string getName() const override;
    [[nodiscard]] std::vector<std::string> getResultHeader() const override;
    [[nodiscard]] std::vector<std::string> getResultFields() const override;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_PARTICIPATIONENTROPY_H
