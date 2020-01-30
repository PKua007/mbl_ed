//
// Created by Piotr Kubala on 21/01/2020.
//

#ifndef MBL_ED_IO_H
#define MBL_ED_IO_H

#endif //MBL_ED_IO_H

#include <string>
#include <vector>

#include "Parameters.h"
#include "analyzer/Analyzer.h"

class IO {
private:
    std::ostream &out;

    void saveOutputToFile(const std::string &header, const std::string &fields, const std::string &outputFilename);

public:
    explicit IO(std::ostream &out) : out{out} { }

    void storeAnalyzerResults(const Parameters &params, const Analyzer &analyzer,
                              const std::vector<std::string> &paramsToPrint, const std::string &fileSignature,
                              const std::string &outputFilename);
    Parameters loadParameters(const std::string &inputFilename, const std::vector<std::string> &overridenParams);
    std::vector<std::string> findEigenenergyFiles(const std::string &directory, const std::string &fileSignature);
};
