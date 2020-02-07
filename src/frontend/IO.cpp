//
// Created by Piotr Kubala on 21/01/2020.
//

#include <filesystem>
#include <iterator>

#include "IO.h"
#include "utils/Utils.h"

void IO::saveOutputToFile(const std::string &header, const std::string &fields, const std::string &outputFilename) {
    bool fileExists = std::ifstream(outputFilename).is_open();
    std::ofstream output(outputFilename, std::fstream::app);
    if (!output)
        die("Cannot open " + outputFilename + " to write analyzer output");
    if (!fileExists)
        output << header << std::endl;

    output << fields << std::endl;
}

std::string IO::stringifyRow(std::vector<std::string> row) const {
    auto csvEscaper = [](std::string &entry) {
        if (entry.find(' ') != std::string::npos)
            entry = "\"" + entry + "\"";
    };

    std::for_each(row.begin(), row.end(), csvEscaper);
    std::ostringstream ostream;
    std::copy(row.begin(), row.end(), std::ostream_iterator<std::string>(ostream, " "));
    return ostream.str();
}

std::string IO::prepareInlineResultHeaderRow(const Analyzer &analyzer,
                                             const std::vector<std::string> &parametersToPrint)
{
    std::vector<std::string> header;
    for (const auto &param : parametersToPrint)
        header.push_back(param);
    auto analyzerHeader = analyzer.getInlineResultsHeader();
    header.insert(header.end(), analyzerHeader.begin(), analyzerHeader.end());
    return this->stringifyRow(header);
}

std::string IO::prepareInlineResultFieldsRow(const Parameters &parameters, const Analyzer &analyzer,
                                             const std::vector<std::string> &parametersToPrint)
{
    std::vector<std::string> fields;
    for (const auto &param : parametersToPrint)
        fields.push_back(parameters.getByName(param));
    auto analyzerFields = analyzer.getInlineResultsFields();
    fields.insert(fields.end(), analyzerFields.begin(), analyzerFields.end());
    return this->stringifyRow(fields);
}


void IO::printInlineAnalyzerResults(const Parameters &params, const Analyzer &analyzer,
                                    const std::vector<std::string> &paramsToPrint)
{
    std::string headerRow = this->prepareInlineResultHeaderRow(analyze
            r, paramsToPrint);
    std::string fieldsRow = this->prepareInlineResultFieldsRow(params, analyzer, paramsToPrint);
    this->logger << std::endl << headerRow << std::endl << fieldsRow << std::endl;
}

void IO::storeAnalyzerResults(const Parameters &params, const Analyzer &analyzer,
                              const std::vector<std::string> &paramsToPrint,
                              std::optional<std::string> inlineResultFilename)
{
    if (inlineResultFilename.has_value()) {
        std::string headerRow = this->prepareInlineResultHeaderRow(analyzer, paramsToPrint);
        std::string fieldsRow = this->prepareInlineResultFieldsRow(params, analyzer, paramsToPrint);
        this->saveOutputToFile(headerRow, fieldsRow, inlineResultFilename.value());
    }

    this->logger << std::endl << "Storing bulk results... " << std::flush;
    analyzer.storeBulkResults(params.getOutputFileSignature());
    this->logger << "done." << std::endl;
}

Parameters IO::loadParameters(const std::string &inputFilename, const std::vector<std::string> &overridenParams) {
    std::ifstream paramsFile(inputFilename);
    if (!paramsFile)
        die("Cannot open " + inputFilename + " to read parameters from.");
    std::stringstream paramsStream;
    paramsStream << paramsFile.rdbuf() << std::endl;
    for (const auto &overridenParam : overridenParams)
        paramsStream << overridenParam << std::endl;

    return Parameters(paramsStream);
}

std::vector<std::string> IO::findEigenenergyFiles(const std::string &directory, const std::string &fileSignature) {
    std::vector<std::string> files;
    for (const auto &entry : std::filesystem::directory_iterator(directory)) {
        const std::filesystem::path &filePath = entry.path();
        std::string prefix = fileSignature + "_";
        if (startsWith(filePath.filename(), prefix) && endsWith(filePath.filename(), "_nrg.bin"))
            files.push_back(filePath);
    }

    return files;
}
