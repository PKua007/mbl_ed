//
// Created by Piotr Kubala on 21/01/2020.
//

#ifndef MBL_ED_IO_H
#define MBL_ED_IO_H

#include <string>
#include <vector>

#include "Parameters.h"
#include "analyzer/Analyzer.h"

/**
 * @brief Some input-output helper methods.
 */
class IO {
private:
    std::ostream &logger;

    void saveOutputToFile(const std::string &header, const std::string &fields, const std::string &outputFilename);
    [[nodiscard]] std::string stringifyRow(std::vector<std::string> row) const;
    std::string prepareInlineResultHeaderRow(const Analyzer &analyzer,
                                             const std::vector<std::string> &parametersToPrint);
    std::string prepareInlineResultFieldsRow(const Parameters &parameters, const Analyzer &analyzer,
                                             const std::vector<std::string> &parametersToPrint);

public:
    /**
     * @param logger Logger used to print some info on storing/restoring stuff.
     */
    explicit IO(std::ostream &logger) : logger{logger} { }

    /**
     * @brief Loads parameters from a file named @a inputFilename in key=value format, however override some params
     * using @a key=value or @a term.key=value strings from @a overridenParams
     */
    Parameters loadParameters(const std::string &inputFilename, const std::vector<std::string> &overridenParams);

    /**
     * @brief Return paths to all files starting with @a fileSignature within @a directory directory
     */
    std::vector<std::string> findEigenenergyFiles(const std::string &directory, const std::string &fileSignature);

    /**
     * @brief Prints inline analyzer result to @a logger from the constructor.
     * @details Format: two rows - header and field values. Last columns are fields from InlineAnalyzerTask -s from
     * @a analyzer, while first few are parameters from @a params specified by @a paramsToPrint
     */
    void printInlineAnalyzerResults(const Parameters &params, const Analyzer &analyzer,
                                    const std::vector<std::string> &paramsToPrint);

    /**
     * @brief Stores BulkAnalyzerTask result to files with prefix given by Parameters::getOutputFileSignature and
     * optionally store inline results to @a inlineResultFilename file.
     * @details Inline results are stores in the same format as in printInlineAnalyzerResults()
     * @param params Parameters for inline results and to generate file signature
     * @param analyzer Analyzer whose results to store
     * @param paramsToPrint which params to include in inline results, see printInlineAnalyzerResults()
     * @param inlineResultFilename if @a std::nullopt no file for inline results will be created, otherwise it will be
     * save to the file with this name
     */
    void storeAnalyzerResults(const Parameters &params, const Analyzer &analyzer,
                              const std::vector<std::string> &paramsToPrint,
                              std::optional<std::string> inlineResultFilename);
};

#endif //MBL_ED_IO_H