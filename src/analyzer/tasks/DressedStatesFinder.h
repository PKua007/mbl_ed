//
// Created by Piotr Kubala on 18/07/2020.
//

#ifndef MBL_ED_DRESSEDSTATESFINDER_H
#define MBL_ED_DRESSEDSTATESFINDER_H


#include "analyzer/BulkAnalyzerTask.h"

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

        friend std::ostream &operator<<(std::ostream &out, const Entry &entry);
    };

    friend std::ostream &operator<<(std::ostream &out, const Entry &entry);

    double relativeMiddleEnergy{};
    double relativeMargin{};
    double coefficientThreshold{};
    std::size_t simulationIdx{};
    std::vector<Entry> result;

public:
    /**
     * @brief Constructs the class.
     * @param relativeMiddleEnergy the middle of the band of normalized (to [0, 1]) eigenenergies
     * @param relativeMargin the width of the band
     * @param coefficientThreshold the threshold, above which one of coefficient of eigenvector has to be to be deemed
     * the dressed state. It has to be higher than sqrt(2)
     */
    DressedStatesFinder(double relativeMiddleEnergy, double relativeMargin, double coefficientThreshold);

    void analyze(const Eigensystem &eigensystem, std::ostream &logger) override;
    [[nodiscard]] std::string getName() const override;

    /**
     * @brief Stores the result, where info about each dressed state found is in one row.
     * @details See Entry 's operator<< to know what is printed in each row.
     */
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_DRESSEDSTATESFINDER_H
