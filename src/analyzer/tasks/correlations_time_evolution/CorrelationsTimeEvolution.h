//
// Created by Piotr Kubala on 17/02/2020.
//

#ifndef MBL_ED_CORRELATIONSTIMEEVOLUTION_H
#define MBL_ED_CORRELATIONSTIMEEVOLUTION_H

#include "analyzer/BulkAnalyzerTask.h"
#include "utils/Assertions.h"
#include "SymmetricMatrix.h"
#include "OccupationEvolution.h"
#include "CorrelationsTimeEntry.h"

class CorrelationsTimeEvolution : public BulkAnalyzerTask {
private:
    struct VectorEvolution {
        FockBase::Vector initialVector;
        std::vector<CorrelationsTimeEntry> timeEntries{};

        [[nodiscard]] std::string getHeader() const;
    };

    std::size_t borderSize{};
    std::vector<VectorEvolution> evolutions{};
    std::vector<double> times{};
    double minTime{};
    double maxTime{};
    std::size_t numSteps{};

    [[nodiscard]] std::size_t getNumberOfSites() const;
    [[nodiscard]] bool hasTimeEntries() const;

public:
    CorrelationsTimeEvolution(double minTime, double maxTime, std::size_t numSteps, std::size_t borderSize,
                              const std::vector<FockBase::Vector> &vectorsToEvolve);

    void analyze(const Eigensystem &eigensystem) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_CORRELATIONSTIMEEVOLUTION_H
