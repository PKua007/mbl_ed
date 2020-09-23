//
// Created by Piotr Kubala on 18/07/2020.
//

#include <sstream>
#include <algorithm>
#include <iterator>

#include "FockVector.h"
#include "utils/Assertions.h"

FockVector::FockVector(const std::string &occupationRepresentation) {
    std::istringstream in(occupationRepresentation);
    std::string token;
    while (std::getline(in, token, '.')) {
        int occupation = std::stoi(token);
        Expects(occupation >= 0);
        this->data.push_back(occupation);
    }
}

std::ostream &operator<<(std::ostream &out, const FockVector &fw) {
    if (fw.empty())
        return out;

    std::copy(fw.data.begin(), fw.data.end() - 1, std::ostream_iterator<int>(out, "."));
    out << fw.data.back();

    return out;
}

FockVector::FockVector(std::size_t sites, const std::string &tag) {
    if (sites == 0)
        return;

    if (tag == "unif") {
        this->data.clear();
        this->data.resize(sites, 1);
    } else if (tag == "dw") {
        if (sites % 2 != 0)
            throw FockVectorParseException("dw Fock vector available only for even number of sites");

        this->data.clear();
        this->data.resize(sites);
        for (std::size_t i{}; i < sites; i += 2)
            this->data[i] = 2;
    } else {
        throw FockVectorParseException("unknown tag: " + tag);
    }
}

FockVector operator+(const FockVector &fw1, const FockVector &fw2) {
    FockVector result;
    result.data.reserve(fw1.size());
    result.data = fw1.data;
    result.data.insert(result.data.end(), fw2.data.begin(), fw2.data.end());
    return result;
}
