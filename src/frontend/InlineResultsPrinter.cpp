//
// Created by Piotr Kubala on 18/01/2020.
//

#include "InlineResultsPrinter.h"

std::string InlineResultsPrinter::stringifyRow(std::vector<std::string> row) const {
    auto csvEscaper = [](std::string &entry) {
        if (entry.find(' ') != std::string::npos)
            entry = "\"" + entry + "\"";
    };

    std::for_each(row.begin(), row.end(), csvEscaper);
    std::ostringstream ostream;
    std::copy(row.begin(), row.end(), std::ostream_iterator<std::string>(ostream, " "));
    return ostream.str();
}

InlineResultsPrinter::InlineResultsPrinter(const Parameters &parameters, const std::vector<std::string> &headerFields,
                                           const std::vector<std::string> &resultsFields,
                                           const std::vector<std::string> &parametersToPrint)
{
    for (const auto &param : parametersToPrint) {
        header.push_back(param);
        fields.push_back(parameters.getByName(param));
    }
    this->fields.insert(this->fields.end(), resultsFields.begin(), resultsFields.end());
    this->header.insert(this->header.end(), headerFields.begin(), headerFields.end());
}
