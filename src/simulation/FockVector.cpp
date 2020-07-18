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
            throw std::runtime_error("dw Fock vector available only for even number of sites");

        this->data.clear();
        this->data.resize(sites);
        for (std::size_t i{}; i < sites; i += 2)
            this->data[i] = 2;
    } else {
        throw std::runtime_error("unknown tag: " + tag);
    }
}
