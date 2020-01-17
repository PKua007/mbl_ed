//
// Created by Piotr Kubala on 28/12/2019.
//

#include <string>
#include <ostream>

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
            throw std::runtime_error("[Parameters::Parameters] Unknown parameter " + key);
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
    out << "number of sites       : " << this->numberOfSites << std::endl;
    out << "number of bosons      : " << this->numberOfBosons << std::endl;
    out << "J                     : " << this->J << std::endl;
    out << "W                     : " << this->W << std::endl;
    out << "U                     : " << this->U << std::endl;
    out << "U1                    : " << this->U1 << std::endl;
    out << "beta                  : " << this->beta << std::endl;
    out << "phi0                  : " << this->phi0 << std::endl;
    out << "usePeriodicBC         : " << (this->usePeriodicBC ? "true" : "false") << std::endl;
    out << "saveEigenenergies     : " << (this->saveEigenenergies ? "true" : "false") << std::endl;
    out << "number of simulations : " << this->numberOfSimulations << std::endl;
    out << "seed                  : " << this->seed << std::endl;
}
