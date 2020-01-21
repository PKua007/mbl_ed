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

Parameters IO::loadParameters(const std::string &inputFilename) {
    std::ifstream input(inputFilename);
    if (!input)
        die("Cannot open " + inputFilename + " to read input parameters");

    Parameters params(input);
    params.print(this->out);
    this->out << std::endl;
    return params;
}

std::vector<std::string> IO::findEigenenergyFiles(const std::string &directory, const std::string &fileSignature) {
    std::vector<std::string> files;
    for (const auto &entry : std::filesystem::directory_iterator(directory)) {
        std::string filename = entry.path();
        std::string prefix = directory + "/" + fileSignature + "_";
        if (startsWith(filename, prefix) && endsWith(filename, "_nrg.dat"))
            files.push_back(filename);
    }

    return files;
}

std::vector<double> IO::loadEigenenergies(const std::string &filename) {
    std::ifstream energiesFile(filename);
    if (!energiesFile)
        die("Cannot open " + filename + " to read eigenenergies from");

    std::vector<double> eigenenergies;
    std::copy(std::istream_iterator<double>(energiesFile), std::istream_iterator<double>(),
              std::back_inserter(eigenenergies));
    return eigenenergies;
}
