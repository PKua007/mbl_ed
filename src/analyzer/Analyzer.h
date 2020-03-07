//
// Created by Piotr Kubala on 16/01/2020.
//

#ifndef MBL_ED_ANALYZER_H
#define MBL_ED_ANALYZER_H


#include <vector>
#include <memory>

#include "simulation/Eigensystem.h"
#include "AnalyzerTask.h"
#include "utils/FileUtils.h"

/**
 * @brief A class which will perform all specified AnalyzerTasks on each Eigensystem passed to be analyzed.
 *
 * The class accepts two types: InlineAnalyzerTask and BulkAnalyzerTask. See each of them for the description.
 */
class Analyzer {
private:
    std::unique_ptr<FileOstreamProvider> ostreamProvider;
    std::vector<std::unique_ptr<AnalyzerTask>> tasks;

public:
    /**
     * @brief Default constructor which uses default instance of FileOstreamProvider (not mocked).
     */
    Analyzer() : Analyzer(std::make_unique<FileOstreamProvider>()) { }

    /**
     * @brief Constructor which accepts user-definde (mocked) FileOstreamProvider for testing.
     */
    explicit Analyzer(std::unique_ptr<FileOstreamProvider> ostreamProvider)
            : ostreamProvider{std::move(ostreamProvider)}
    { }

    /**
     * @brief Appends another task to the Analyzer.
     */
    void addTask(std::unique_ptr<AnalyzerTask> task);

    /**
     * @brief Performs all analzyer tasks added by Analyzer::addTask on this @a eigensystem.
     */
    void analyze(const Eigensystem &eigensystem, std::ostream &logger);

    /**
     * @brief Returns a vector of names of fields imploded from all InlineAnalyzerTask -s. The order is the same
     * as in Analyzer::getInlineResultsFields.
     */
    [[nodiscard]] std::vector<std::string> getInlineResultsHeader() const;

    /**
     * @brief Returns a vector of values of fields imploded from all InlineAnalyzerTask -s. The order is the same
     * as in Analyzer::getInlineResultsHeader.
     */
    [[nodiscard]] std::vector<std::string> getInlineResultsFields() const;

    /**
     * @brief Saves results of BulkAnalyzerTask-s to files.
     *
     * The name of files are given by `fileSignature_name.txt`, where `name` is fetched using AnalyzerTask::getName.
     * Note, that using Analyzer::Analyzer(std::unique_ptr<FileOstreamProvider>) one can avoid actually creating files.
     * @param fileSignature
     */
    void storeBulkResults(const std::string &fileSignature) const;
};


#endif //MBL_ED_ANALYZER_H
