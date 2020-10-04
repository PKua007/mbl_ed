//
// Created by Piotr Kubala on 02/10/2020.
//

#ifndef MBL_ED_EIGENSTATEOBSERVABLES_H
#define MBL_ED_EIGENSTATEOBSERVABLES_H

#include <vector>
#include <string>
#include <memory>

#include "analyzer/BulkAnalyzerTask.h"
#include "core/PrimaryObservable.h"
#include "core/SecondaryObservable.h"
#include "utils/FileUtils.h"

/**
 * @brief A BulkAnalyzerTask calculating values of given observables for epsilons in a specified number of bins.
 * @details For each observable value, we calculate and average over the epsilons landing in a given bin for a single
 * eigensystem, and then this average is a single "experiment", which is then averaged over multiple eigensystems
 * and its error is calculated
 */
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

    std::size_t eigensystemIdx{};
    std::string fileSignature;
    std::unique_ptr<FileOstreamProvider> ostreamProvider;

    static auto calculateBinRange(size_t binIdx, size_t numBins);

    [[nodiscard]] std::string generateStoredObservablesHeader() const;
    [[nodiscard]] std::size_t countStoredObservableValues() const;
    [[nodiscard]] arma::mat calculateObservables(const Eigensystem &eigensystem) const;

    [[nodiscard]] std::vector<double>
    calculateMeanObservableValuesForBand(const arma::mat &observables,
                                         const std::vector<std::size_t> &bandIndices) const;

public:
    EigenstateObservables(std::size_t numOfBins, std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables,
                          std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables,
                          std::vector<std::shared_ptr<Observable>> storedObservables);

    void startStoringObservables(std::string fileSignature_, std::unique_ptr<FileOstreamProvider> ostreamProvider_
                                 = std::make_unique<FileOstreamProvider>())
    {
        this->ostreamProvider = std::move(ostreamProvider_);
        this->fileSignature = std::move(fileSignature_);
    };

    void analyze(const Eigensystem &eigensystem, Logger &logger) override;
    [[nodiscard]] std::string getName() const override;

    /**
     * @brief Each line in out is the entry for subsequent bins with format: [bin start] [observable 1] [obs 1 error]
     * [obs 2] [obs 2 error], ...
     * @details Moreover, first line is a header: binStart [obs 1 name] d[obs1 name] [obs 2 name] d[obs 2 name] ...
     */
    void storeResult(std::ostream &out) const override;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_EIGENSTATEOBSERVABLES_H
