//
// Created by Piotr Kubala on 29/01/2020.
//

#ifndef MBL_ED_INVERSEPARTICIPATIONRATIO_H
#define MBL_ED_INVERSEPARTICIPATIONRATIO_H

#include "analyzer/BulkAnalyzerTask.h"
#include "analyzer/BandExtractor.h"

/**
 * @brief BulkAnalyzerTask which plots inverse participation ratio for indivudial eigenvectors vs corresponding
 * eigenenergies.
 *
 * For a given eigenvector \f$v^i\f$, ipr is defined as \f$\sum_{i=0}^\text{dim}(v^i)^4\f$.
 */
class InverseParticipationRatio : public BulkAnalyzerTask {
private:
    /**
     * @brief X, Y struct for each point on the plot
     */
    struct Entry {
        Entry() = default;
        Entry(double energy, double ipr) : energy(energy), ipr(ipr) { }

        double energy{};
        double ipr{};

        friend std::ostream &operator<<(std::ostream &out, const Entry &entry);
    };

    friend std::ostream &operator<<(std::ostream &out, const Entry &entry);

    BandExtractor extractor;
    std::vector<Entry> entries{};

public:
    /**
     * @brief The constructor, which specifies the band of eigenenergies normalized to [0, 1] from which we will take
     * points to the plot
     * @param relativeMiddleEnergy the middle point of the band
     * @param relativeMargin the width of the band
     */
    explicit InverseParticipationRatio(BandExtractor::Range range) : extractor(std::move(range), "Ipr ") { };

    /**
     * @brief Adds points to individual eigenvectors vs corresponding eigenenergies plot from given @a eigensystem.
     *
     * Invoking multiple times will just append more points.
     * @param eigensystem Eigensystem to add points from
     * @param logger unused
     */
    void analyze(const Eigensystem &eigensystem, Logger &logger) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_INVERSEPARTICIPATIONRATIO_H
