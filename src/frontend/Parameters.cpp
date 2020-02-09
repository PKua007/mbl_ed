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
        if (key == "lattice.numberOfSites")
            this->numberOfSites = config.getUnsignedLong("lattice.numberOfSites");
        else if (key == "lattice.numberOfBosons")
            this->numberOfBosons = config.getUnsignedLong("lattice.numberOfBosons");
        else if (key == "lattice.usePeriodicBC")
            this->usePeriodicBC = config.getBoolean("lattice.usePeriodicBC");

        else if (key == "hamiltonian.J")
            this->J = config.getDouble("hamiltonian.J");
        else if (key == "hamiltonian.W")
            this->W = config.getDouble("hamiltonian.W");
        else if (key == "hamiltonian.U")
            this->U = config.getDouble("hamiltonian.U");
        else if (key == "hamiltonian.U1")
            this->U1 = config.getDouble("hamiltonian.U1");
        else if (key == "hamiltonian.beta")
            this->beta = config.getDouble("hamiltonian.beta");
        else if (key == "hamiltonian.phi0")
            this->phi0 = config.getDouble("hamiltonian.phi0");

        else if (key == "simulation.averagingModel")
            this->averagingModel = config.getString("simulation.averagingModel");
        else if (key == "simulation.calculateEigenvectors")
            this->calculateEigenvectors = config.getBoolean("simulation.calculateEigenvectors");
        else if (key == "simulation.saveEigenenergies")
            this->saveEigenenergies = config.getBoolean("simulation.saveEigenenergies");
        else if (key == "simulation.from")
            this->from = config.getUnsignedLong("simulation.from");
        else if (key == "simulation.to")
            this->to = config.getUnsignedLong("simulation.to");
        else if (key == "simulation.totalSimulations")
            this->totalSimulations = config.getUnsignedLong("simulation.totalSimulations");
        else if (key == "simulation.seed")
            this->seed = config.getUnsignedLong("simulation.seed");

        else
            throw UnknownParameterException("Unknown parameter " + key);
    }

    this->autocompleteAndValidate();
}

void Parameters::autocompleteAndValidate() {
    // Empty totalSimulations is taken from to and vice verse, but at least one of them must be specified
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
    out << "averagingModel        : " << this->averagingModel << std::endl;
    out << "usePeriodicBC         : " << (this->usePeriodicBC ? "true" : "false") << std::endl;
    out << "calculateEigenvectors : " << (this->calculateEigenvectors ? "true" : "false") << std::endl;
    out << "saveEigenenergies     : " << (this->saveEigenenergies ? "true" : "false") << std::endl;
    out << "from                  : " << this->from << std::endl;
    out << "to                    : " << this->to << std::endl;
    out << "totalSimulations      : " << this->totalSimulations << std::endl;
    out << "seed                  : " << this->seed << std::endl;
    out << std::endl;
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
        return this->doubleToString(this->phi0);
    else if (name == "averagingModel")
        return this->averagingModel;
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