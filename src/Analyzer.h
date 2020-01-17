//
// Created by Piotr Kubala on 16/01/2020.
//

#ifndef MBL_ED_ANALYZER_H
#define MBL_ED_ANALYZER_H


#include <vector>
#include <memory>

#include "AnalyzerTask.h"
#include "FileUtils.h"

class Analyzer {
private:
    std::unique_ptr<FileOstreamProvider> ostreamProvider;
    std::vector<std::unique_ptr<AnalyzerTask>> tasks;

public:
    Analyzer() : Analyzer(std::make_unique<FileOstreamProvider>()) { }
    explicit Analyzer(std::unique_ptr<FileOstreamProvider> ostreamProvider)
            : ostreamProvider{std::move(ostreamProvider)}
    { }

    void addTask(std::unique_ptr<AnalyzerTask> task);
    void analyze(const std::vector<double> &eigenenergies);
    [[nodiscard]] std::vector<std::string> getInlineResultsHeader() const;
    [[nodiscard]] std::vector<std::string> getInlineResultsFields() const;
    void storeBulkResults(const std::string &fileSignature) const;
};


#endif //MBL_ED_ANALYZER_H
