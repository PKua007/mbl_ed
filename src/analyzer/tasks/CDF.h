//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_CDF_H
#define MBL_ED_CDF_H


#include "analyzer/BulkAnalyzerTask.h"

class CDF : public BulkAnalyzerTask {
private:
    std::vector<std::vector<double>> cdfTable;

    [[nodiscard]] std::vector<double> getNormalizedEigenenergies(const std::vector<double> &eigenenergies) const;

public:
    explicit CDF(std::size_t bins);

    void analyze(const Eigensystem &eigensystem) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_CDF_H
