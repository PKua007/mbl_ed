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

    private:
        std::string getInitialVectorSignature() const;
    };

    std::size_t marginSize{};
    std::vector<VectorEvolution> vectorEvolutions{};
    double minTime{};
    double maxTime{};
    std::size_t numSteps{};

    [[nodiscard]] std::size_t getNumberOfSites() const;
    [[nodiscard]] bool numberOfSitesDetermined() const;
    void prepareTimeEntriesForNumberOfSites(std::size_t numberOfSites);

public:
    CorrelationsTimeEvolution(double minTime, double maxTime, std::size_t numSteps, std::size_t marginSize,
                              const std::vector<FockBase::Vector> &vectorsToEvolve);

    void analyze(const Eigensystem &eigensystem) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_CORRELATIONSTIMEEVOLUTION_H
