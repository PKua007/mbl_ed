//
// Created by Piotr Kubala on 31/07/2020.
//

#ifndef MBL_ED_BULKMEANGAPRATIO_H
#define MBL_ED_BULKMEANGAPRATIO_H


#include "analyzer/BulkAnalyzerTask.h"

#include "utils/Assertions.h"

class BulkMeanGapRatio : public BulkAnalyzerTask {
private:
    std::vector<std::vector<double>> gapRatios;

public:
    explicit BulkMeanGapRatio(std::size_t numBins) : gapRatios(numBins) { Expects(numBins > 0); }

    void analyze(const Eigensystem &eigensystem, std::ostream &logger) override;
    [[nodiscard]] std::string getName() const override { return "mgrs"; }
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_BULKMEANGAPRATIO_H
