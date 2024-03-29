//
// Created by Piotr Kubala on 18/07/2020.
//

#ifndef MBL_ED_DRESSEDSTATESFINDER_H
#define MBL_ED_DRESSEDSTATESFINDER_H


#include "analyzer/BulkAnalyzerTask.h"
#include "analyzer/BandExtractor.h"

/**
 * @brief Analyzer class, which finds dressed states in a given energy band.
 * @details Namely, it looks for eigenvectors, where one of the coefficients in Fock basis is above a given threshold.
 */
class DressedStatesFinder : public BulkAnalyzerTask {
private:
    /**
     * @brief Helper class storing the entry in each row of the output.
     */
    struct Entry {
        std::size_t simulationIdx{};
        std::string vector;
        double energy{};
        double coefficient{};

        void storeState(std::ostream &binaryOut) const;
        void restoreState(std::istream &binaryIn);

        friend std::ostream &operator<<(std::ostream &out, const Entry &entry);
    };

    friend std::ostream &operator<<(std::ostream &out, const Entry &entry);

    BandExtractor extractor;
    double coefficientThreshold{};
    std::size_t simulationIdx{};
    std::vector<Entry> result;

public:
    /**
     * @brief Constructs the class.
     * @param range the range to choose eigenstates from
     * @param coefficientThreshold the threshold, above which one of coefficient of eigenvector has to be to be deemed
     * the dressed state. It has to be higher than sqrt(2)
     */
    DressedStatesFinder(double coefficientThreshold, BandExtractor::Range range);

    void analyze(const Eigensystem &eigensystem, Logger &logger) override;
    [[nodiscard]] std::string getName() const override;

    /**
     * @brief Stores the result, where info about each dressed state found is in one row.
     * @details See Entry 's operator<< to know what is printed in each row.
     */
    void storeResult(std::ostream &out) const override;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_DRESSEDSTATESFINDER_H
