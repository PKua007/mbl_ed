//
// Created by Piotr Kubala on 21/01/2020.
//

#include <filesystem>
#include <iterator>

#include "IO.h"
#include "utils/Utils.h"

/**
 * @brief Saves inline results' @a header + @a fields to @a outputFilename. The header is printed only if file didn't
 * exist.
 */
void IO::saveOutputToFile(const std::string &header, const std::string &fields, const std::string &outputFilename) {
    bool fileExists = std::ifstream(outputFilename).is_open();
    std::ofstream output(outputFilename, std::fstream::app);
    if (!output)
        die("Cannot open " + outputFilename + " to write analyzer output");
    if (!fileExists)
        output << header << std::endl;

    output << fields << std::endl;
}

/**
 * @brief Implodes a vector of strings into a row with spaces escaped by " "
 */
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

/**
 * @brief Takes @a analyzer and prepares the whole header for inline results.
 * @param parametersToPrint parameters from Parameters to also include (the first columns)
 */
std::string IO::prepareInlineResultHeaderRow(const std::vector<std::string> &parametersToPrint,
        const std::vector<std::string> &additionalFields)
{
    std::vector<std::string> header;
    for (const auto &param : parametersToPrint)
        header.push_back(param);
    header.insert(header.end(), additionalFields.begin(), additionalFields.end());
    return this->stringifyRow(header);
}

/**
 * @brief Takes @a analyzer and prepares the whole fields row for inline results.
 * @param parametersToPrint parameters from Parameters to also include (the first columns)
 */
std::string IO::prepareInlineResultFieldsRow(const Parameters &parameters,
        const std::vector<std::string> &parametersToPrint, const std::vector<std::string> &additionalFields)
{
    std::vector<std::string> fields;
    for (const auto &param : parametersToPrint)
        fields.push_back(parameters.getByName(param));
    fields.insert(fields.end(), additionalFields.begin(), additionalFields.end());
    return this->stringifyRow(fields);
}


void IO::printInlineResults(const Parameters &params, const std::vector<std::string> &paramsToPrint,
        const std::vector<std::string> &additionalHeader, const std::vector<std::string> &additionalFields)
{
    std::string headerRow = this->prepareInlineResultHeaderRow(paramsToPrint, additionalHeader);
    std::string fieldsRow = this->prepareInlineResultFieldsRow(params, paramsToPrint, additionalFields);
    this->logger << std::endl << headerRow << std::endl << fieldsRow << std::endl;
}

void IO::storeAnalyzerResults(const Parameters &params, const Analyzer &analyzer,
                              const std::vector<std::string> &paramsToPrint,
                              std::optional<std::string> inlineResultFilename)
{
    if (inlineResultFilename.has_value()) {
        this->storeInlineResults(params, paramsToPrint, analyzer.getInlineResultsHeader(),
                                 analyzer.getInlineResultsFields(), *inlineResultFilename);
    }

    this->logger << std::endl << "Storing bulk results... " << std::flush;
    analyzer.storeBulkResults(params.getOutputFileSignatureWithRange());
    this->logger << "done." << std::endl;
}

Parameters IO::loadParameters(const std::string &inputFilename, const std::vector<std::string> &overridenParams) {
    std::ifstream paramsFile(inputFilename);
    if (!paramsFile)
        die("Cannot open " + inputFilename + " to read parameters from.");
    std::stringstream paramsStream;
    paramsStream << paramsFile.rdbuf() << std::endl;

    for (const auto &overridenParam : overridenParams) {
        std::size_t dotPos = overridenParam.find('.');
        std::size_t equalPos = overridenParam.find('=');
        if (equalPos == std::string::npos)
            die("Malformed overriden param. Use: [param name]=[value] or [term name].[param name]=[value]");

        if (dotPos == std::string::npos || dotPos > equalPos) {
            // Overriding general param
            // 1) no '.', or
            // 2) dot after '=', meaning that it is a part of the value, eg. someParam=1.5
            paramsStream << "[general]" << std::endl << overridenParam << std::endl;
        } else {
            // Overriding hamiltonian term param
            if (dotPos == 0 || dotPos == overridenParam.size() - 1)
                die("Malformed overriden param. Use: [param name]=[value] or [term name].[param name]=[value]");
            std::string termName = overridenParam.substr(0, dotPos);
            std::string paramName = overridenParam.substr(dotPos + 1);
            paramsStream << "[term." << termName << "]" << std::endl << paramName << std::endl;
        }
    }

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

void IO::storeInlineResults(const Parameters &parameters, const std::vector<std::string> &paramsToStore,
                            const std::vector<std::string> &header, const std::vector<std::string> &fields,
                            const std::string &outputFilename) {
    std::string headerRow = this->prepareInlineResultHeaderRow(paramsToStore, header);
    std::string fieldsRow = this->prepareInlineResultFieldsRow(parameters, paramsToStore, fields);
    this->saveOutputToFile(headerRow, fieldsRow, outputFilename);
}
