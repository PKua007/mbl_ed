//
// Created by Piotr Kubala on 17/02/2020.
//

#ifndef MBL_ED_CORRELATIONSTIMEEVOLUTION_H
#define MBL_ED_CORRELATIONSTIMEEVOLUTION_H

#include "analyzer/BulkAnalyzerTask.h"
#include "utils/Assertions.h"
#include "SymmetricMatrix.h"
#include "OccupationEvolution.h"

class CorrelationsTimeEvolution : public BulkAnalyzerTask {
public:
    enum TimeScaleType {
        Linear,
        Logarithmic
    };

private:
    struct Correlations {
        std::size_t distance{};
        double G{};

        void addObservables(const OccupationEvolution::Observables &observables, std::size_t borderSize);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string getValue(std::size_t meanEntries) const;
    };

    struct OnsiteFluctuations {
        std::size_t i{};
        double rho{};

        void addObservables(const OccupationEvolution::Observables &observables);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string getValue(std::size_t meanEntries) const;
    };

    struct TimeEntry {
    private:
        double t{};
        double x{};
        std::vector<Correlations> correlations{};
        std::vector<Correlations> borderlessCorrelations{};
        std::vector<OnsiteFluctuations> onsiteFluctuations{};

        std::size_t borderSize{};
        std::size_t meanEntries{};

    public:
        TimeEntry() = default;
        TimeEntry(double t, std::size_t borderSize, std::size_t numberOfSites);

        [[nodiscard]] std::size_t getNumberOfSites() const;
        void addObservables(const OccupationEvolution::Observables &observables);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string getValue() const;
    };

    struct VectorEvolution {
        FockBase::Vector initialVector;
        std::vector<TimeEntry> timeEntries{};

        [[nodiscard]] std::string getHeader() const;
    };

    std::size_t borderSize{};
    std::vector<VectorEvolution> evolutions{};
    std::vector<double> times{};
    std::size_t meanEntries{};

    [[nodiscard]] std::size_t getNumberOfSites() const;
    [[nodiscard]] bool hasTimeEntries() const;

public:
    CorrelationsTimeEvolution(double minTime, double maxTime, std::size_t numSteps, TimeScaleType timeScaleType,
                              std::size_t borderSize, const std::vector<FockBase::Vector> &vectorsToEvolve);

    void analyze(const Eigensystem &eigensystem) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_CORRELATIONSTIMEEVOLUTION_H
