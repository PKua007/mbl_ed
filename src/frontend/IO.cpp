//
// Created by Piotr Kubala on 21/01/2020.
//

#include <filesystem>

#include "IO.h"
#include "utils/Utils.h"
#include "InlineResultsPrinter.h"

void IO::saveOutputToFile(const std::string &header, const std::string &fields, const std::string &outputFilename) {
    bool fileExists = std::ifstream(outputFilename).is_open();
    std::ofstream output(outputFilename, std::fstream::app);
    if (!output)
        die("Cannot open " + outputFilename + " to write analyzer output");
    if (!fileExists)
        output << header << std::endl;

    output << fields << std::endl;
}

void IO::storeAnalyzerResults(const Parameters &params, const Analyzer &analyzer,
                              const std::vector<std::string> &paramsToPrint, const std::string &fileSignature,
                              const std::string &outputFilename)
{
    InlineResultsPrinter resultsPrinter(params, analyzer.getInlineResultsHeader(),
                                        analyzer.getInlineResultsFields(), paramsToPrint);
    this->out << std::endl << resultsPrinter.getHeader() << std::endl << resultsPrinter.getFields() << std::endl;

    if (!outputFilename.empty())
        saveOutputToFile(resultsPrinter.getHeader(), resultsPrinter.getFields(), outputFilename);

    this->out << std::endl << "Storing bulk results... " << std::flush;
    analyzer.storeBulkResults(fileSignature);
    this->out << "done." << std::endl;
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
