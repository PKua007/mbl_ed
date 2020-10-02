//
// Created by Piotr Kubala on 02/10/2020.
//

#ifndef MBL_ED_EIGENSTATEOBSERVABLES_H
#define MBL_ED_EIGENSTATEOBSERVABLES_H

#include <vector>
#include <string>
#include <memory>

#include "analyzer/BulkAnalyzerTask.h"
#include "evolution/PrimaryObservable.h"
#include "evolution/SecondaryObservable.h"

class EigenstateObservables : public BulkAnalyzerTask {
private:
    struct BinEntry : public Restorable {
        std::vector<std::vector<double>> observableValues;

        void storeState(std::ostream &binaryOut) const override;

        void joinRestoredState(std::istream &binaryIn) override;

        void clear() override;
    };

    std::string header;
    std::vector<BinEntry> binEntries;
    std::size_t numValues{};

    std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables;
    std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables;
    std::vector<std::shared_ptr<Observable>> storedObservables;

    [[nodiscard]] std::string generateStoredObservablesHeader() const;
    [[nodiscard]] std::size_t countStoredObservableValues() const;

public:
    EigenstateObservables(std::size_t numOfBins, std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables,
                          std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables,
                          std::vector<std::shared_ptr<Observable>> storedObservables);

    void analyze(const Eigensystem &eigensystem, Logger &logger) override;
    std::string getName() const override;
    void storeResult(std::ostream &out) const override;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_EIGENSTATEOBSERVABLES_H
