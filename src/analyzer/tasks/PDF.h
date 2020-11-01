//
// Created by Piotr Kubala on 01/11/2020.
//

#ifndef MBL_ED_PDF_H
#define MBL_ED_PDF_H

#include "analyzer/BulkAnalyzerTask.h"

class PDF : public BulkAnalyzerTask {
private:
    std::vector<std::vector<double>> pdfTable;

public:
    /**
     * @brief Constructs the class.
     * @param bins Number of bins that will be filled.
     */
    explicit PDF(std::size_t bins);

    /**
     * @brief Takes normalized eigenenergies from the Eigensystem and inserts it into the PDF histogram.
     *
     * Averaging multiple eigensystems: it produces a separate, normalized graph for each eigensystem and then
     * averages all this separate graphs (note, that the final graph is also normalized).
     */
    void analyze(const Eigensystem &binValue, Logger &logger) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_PDF_H
