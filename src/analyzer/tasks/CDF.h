//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_CDF_H
#define MBL_ED_CDF_H


#include "analyzer/BulkAnalyzerTask.h"

/**
 * @brief BulkAnalyzerTask, which draws cumulative distribution function for normalized eigenenergies (from [0, 1])
 */
class CDF : public BulkAnalyzerTask {
private:
    std::vector<std::vector<double>> cdfTable;

public:
    /**
     * @brief Constructs the class.
     * @param bins Number of bins that will be filled.
     */
    explicit CDF(std::size_t bins);

    /**
     * @brief Takes normalized eigenenergies from the Eigensystem and inserts it into the histogram.
     *
     * Averaging multiple eigensystems: it produces a separate, normalized graph for each eigensystem and then
     * averages all this separate graphs (note, that the final graph is also normalized).
     * @param eigensystem Eigensystem to draw histogram from
     */
    void analyze(const Eigensystem &eigensystem) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_CDF_H
