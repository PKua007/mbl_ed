//
// Created by Piotr Kubala on 18/01/2020.
//

#ifndef MBL_ED_INLINERESULTSPRINTER_H
#define MBL_ED_INLINERESULTSPRINTER_H


#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <iterator>

#include "Parameters.h"

class InlineResultsPrinter {
private:
    std::vector<std::string> header;
    std::vector<std::string> fields;

    [[nodiscard]] std::string stringifyRow(std::vector<std::string> row) const;

public:
    InlineResultsPrinter(const Parameters &parameters, const std::vector<std::string> &headerFields,
                         const std::vector<std::string> &resultsFields,
                         const std::vector<std::string> &parametersToPrint);

    [[nodiscard]] std::string getHeader() const { return this->stringifyRow(this->header); }
    [[nodiscard]] std::string getFields() const { return this->stringifyRow(this->fields); }
};


#endif //MBL_ED_INLINERESULTSPRINTER_H
