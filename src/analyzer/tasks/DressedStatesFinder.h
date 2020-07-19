//
// Created by Piotr Kubala on 18/07/2020.
//

#ifndef MBL_ED_DRESSEDSTATESFINDER_H
#define MBL_ED_DRESSEDSTATESFINDER_H


#include "analyzer/BulkAnalyzerTask.h"

class DressedStatesFinder : public BulkAnalyzerTask {
private:
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
    DressedStatesFinder(double relativeMiddleenergy, double relativeMargin, double coefficientThreshold);

    void analyze(const Eigensystem &eigensystem, std::ostream &logger) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_DRESSEDSTATESFINDER_H
