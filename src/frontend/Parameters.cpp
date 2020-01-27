//
// Created by Piotr Kubala on 28/12/2019.
//

#include <string>
#include <ostream>
#include <sstream>

#include "Parameters.h"
#include "utils/Config.h"
#include "utils/Assertions.h"

Parameters::Parameters(std::istream &input) {
    auto config = Config::parse(input, '=', true);
    for (const auto &key : config.getKeys()) {
        if (key == "numberOfSites")
            this->numberOfSites = config.getUnsignedLong("numberOfSites");
        else if (key == "numberOfBosons")
            this->numberOfBosons = config.getUnsignedLong("numberOfBosons");
        else if (key == "J")
            this->J = config.getDouble("J");
        else if (key == "W")
            this->W = config.getDouble("W");
        else if (key == "U")
            this->U = config.getDouble("U");
        else if (key == "U1")
            this->U1 = config.getDouble("U1");
        else if (key == "beta")
            this->beta = config.getDouble("beta");
        else if (key == "phi0")
            this->phi0 = config.getString("phi0");
        else if (key == "usePeriodicBC")
            this->usePeriodicBC = config.getBoolean("usePeriodicBC");
        else if (key == "calculateEigenvectors")
            this->calculateEigenvectors = config.getBoolean("calculateEigenvectors");
        else if (key == "saveEigenenergies")
            this->saveEigenenergies = config.getBoolean("saveEigenenergies");
        else if (key == "from")
            this->from = config.getUnsignedLong("from");
        else if (key == "to")
            this->to = config.getUnsignedLong("to");
        else if (key == "totalSimulations")
            this->totalSimulations = config.getUnsignedLong("totalSimulations");
        else if (key == "seed")
            this->seed = config.getUnsignedLong("seed");
        else
            throw UnknownParameterException("Unknown parameter " + key);
    }

    this->autocompleteAndValidate();
}

void Parameters::autocompleteAndValidate() {
    Validate(this->totalSimulations > 0 || this->to > 0);

    if (this->to == 0)
        this->to = this->totalSimulations;
    else if (this->totalSimulations == 0)
        this->totalSimulations = this->to;

    Validate(this->from < this->to);
    Validate(this->to <= this->totalSimulations);
    Validate(numberOfSites > 0);
    Validate(numberOfBosons > 0);
    Validate(J >= 0);
    Validate(W >= 0);
    Validate(U >= 0);
    Validate(U1 >= 0);
}

void Parameters::print(std::ostream &out) const {
    out << "numberOfSites         : " << this->numberOfSites << std::endl;
    out << "numbeOfBosons         : " << this->numberOfBosons << std::endl;
    out << "J                     : " << this->J << std::endl;
    out << "W                     : " << this->W << std::endl;
    out << "U                     : " << this->U << std::endl;
    out << "U1                    : " << this->U1 << std::endl;
    out << "beta                  : " << this->beta << std::endl;
    out << "phi0                  : " << this->phi0 << std::endl;
    out << "usePeriodicBC         : " << (this->usePeriodicBC ? "true" : "false") << std::endl;
    out << "calculateEigenvectors : " << (this->calculateEigenvectors ? "true" : "false") << std::endl;
    out << "saveEigenenergies     : " << (this->saveEigenenergies ? "true" : "false") << std::endl;
    out << "from                  : " << this->from << std::endl;
    out << "to                    : " << this->to << std::endl;
    out << "totalSimulations      : " << this->totalSimulations << std::endl;
    out << "seed                  : " << this->seed << std::endl;
}

std::string Parameters::getByName(const std::string &name) const {
    if (name == "numberOfSites")
        return this->doubleToString(this->numberOfSites);
    else if (name == "numberOfBosons")
        return this->doubleToString(this->numberOfBosons);
    else if (name == "J")
        return this->doubleToString(this->J);
    else if (name == "W")
        return this->doubleToString(this->W);
    else if (name == "U")
        return this->doubleToString(this->U);
    else if (name == "U1")
        return this->doubleToString(this->U1);
    else if (name == "beta")
        return this->doubleToString(this->beta);
    else if (name == "phi0")
        return this->phi0;
    else if (name == "usePeriodicBC")
        return this->usePeriodicBC ? "true" : "false";
    else if (name == "calculateEigenvectors")
        return this->calculateEigenvectors ? "true" : "false";
    else if (name == "saveEigenenergies")
        return this->saveEigenenergies ? "true" : "false";
    else if (name == "from")
        return std::to_string(this->from);
    else if (name == "to")
        return std::to_string(this->to);
    else if (name == "totalSimulations")
        return std::to_string(this->totalSimulations);
    else if (name == "seed")
        return std::to_string(this->seed);
    else
        throw UnknownParameterException("Unknown parameter " + name);
}

std::string Parameters::doubleToString(double d) const {
    std::ostringstream ostream;
    ostream << d;
    return ostream.str();
}

std::string Parameters::getOutputFileSignature() const {
    std::ostringstream filename;
    filename << "J." << this->J << "_U." << this->U << "_U1." << this->U1 << "_N." << this->numberOfBosons;
    filename << "_K." << this->numberOfSites  << "_beta." << this->beta << "_phi0." << this->phi0;
    return filename.str();
}