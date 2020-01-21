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
        else if (key == "saveEigenenergies")
            this->saveEigenenergies = config.getBoolean("saveEigenenergies");
        else if (key == "numberOfSimulations")
            this->numberOfSimulations = config.getUnsignedLong("numberOfSimulations");
        else if (key == "seed")
            this->seed = config.getUnsignedLong("seed");
        else
            throw UnknownParameterException("Unknown parameter " + key);
    }

    this->validate();
}

void Parameters::validate() const {
    Validate(numberOfSites > 0);
    Validate(numberOfBosons > 0);
    Validate(J >= 0);
    Validate(W >= 0);
    Validate(U >= 0);
    Validate(U1 >= 0);
    Validate(numberOfSimulations > 0);
}

void Parameters::print(std::ostream &out) const {
    out << "numberOfSites       : " << this->numberOfSites << std::endl;
    out << "numbeOfBosons       : " << this->numberOfBosons << std::endl;
    out << "J                   : " << this->J << std::endl;
    out << "W                   : " << this->W << std::endl;
    out << "U                   : " << this->U << std::endl;
    out << "U1                  : " << this->U1 << std::endl;
    out << "beta                : " << this->beta << std::endl;
    out << "phi0                : " << this->phi0 << std::endl;
    out << "usePeriodicBC       : " << (this->usePeriodicBC ? "true" : "false") << std::endl;
    out << "saveEigenenergies   : " << (this->saveEigenenergies ? "true" : "false") << std::endl;
    out << "numberOfSimulations : " << this->numberOfSimulations << std::endl;
    out << "seed                : " << this->seed << std::endl;
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
    else if (name == "saveEigenenergies")
        return this->saveEigenenergies ? "true" : "false";
    else if (name == "numberOfSimulations")
        return std::to_string(this->numberOfSimulations);
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