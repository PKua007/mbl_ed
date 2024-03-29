//
// Created by Piotr Kubala on 28/12/2019.
//

#include <string>
#include <ostream>
#include <sstream>

#include "Parameters.h"
#include "utils/Config.h"
#include "utils/Assertions.h"
#include "utils/Utils.h"

Parameters::Parameters(std::istream &input) {
    // First we are looking for parameters from [general sections]
    auto config = Config::parse(input, '=', true);
    if (!config.hasRootSection("general"))
        throw ParametersParseException("Found no general parameters.");
    
    auto generalConfig = config.fetchSubconfig("general");
    for (const auto &key : generalConfig.getKeys()) {
        if (key == "N")
            this->N = generalConfig.getUnsignedLong("N");
        else if (key == "K")
            this->K = generalConfig.getUnsignedLong("K");
        else if (key == "usePeriodicBC")
            this->usePeriodicBC = generalConfig.getBoolean("usePeriodicBC");
        else if (key == "averagingModel")
            this->averagingModel = generalConfig.getString("averagingModel");
        else if (key == "calculateEigenvectors")
            this->calculateEigenvectors = generalConfig.getBoolean("calculateEigenvectors");
        else if (key == "saveEigenenergies")
            this->saveEigenenergies = generalConfig.getBoolean("saveEigenenergies");
        else if (key == "saveEigenstates")
            this->saveEigenstates = generalConfig.getBoolean("saveEigenstates");
        else if (key == "storeFormat")
            this->storeFormat = generalConfig.getString("storeFormat");
        else if (key == "from")
            this->from = generalConfig.getUnsignedLong("from");
        else if (key == "to")
            this->to = generalConfig.getUnsignedLong("to");
        else if (key == "totalSimulations")
            this->totalSimulations = generalConfig.getUnsignedLong("totalSimulations");
        else if (key == "seed")
            this->seed = generalConfig.getUnsignedLong("seed");
        else if (key == "splitWorkload")
            this->splitWorkload = generalConfig.getBoolean("splitWorkload");
        else if (key == "secureSimulationState")
            this->secureSimulationState = generalConfig.getBoolean("secureSimulationState");
        else
            throw ParametersParseException("Unknown parameter " + key);
    }

    // Now we are parsing all hamiltonian terms - [term.termName] sections
    if (!config.hasRootSection("term"))
        throw ParametersParseException("No hamiltonian terms were specified!");

    Config terms = config.fetchSubconfig("term");
    if (terms.hasRootSection(""))
        throw ParametersParseException("Found naked section [term]. Expecting only [term.nameOfTerm] sections.");

    for (const auto &term : terms.getRootSections())
        this->hamiltonianTerms[term] = terms.fetchSubconfig(term);

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
    Validate(N > 0);
    Validate(K > 0);
    if (this->splitWorkload)
        ValidateMsg(this->secureSimulationState, "If splitWorkoload = true, secureSimulationState should also be true");

    ValidateMsg(!this->saveEigenstates || this->saveEigenenergies,
                "Eigenstates cannot be stored without eigenenergies");
    ValidateMsg(!this->saveEigenstates || this->calculateEigenvectors, "Eigenvectors must be calculated to be stored");
}

void Parameters::printGeneral(std::ostream &out) const {
    out << "General parameters:" << std::endl;
    out << "N                     : " << this->N << std::endl;
    out << "K                     : " << this->K << std::endl;
    out << "averagingModel        : " << this->averagingModel << std::endl;
    out << "usePeriodicBC         : " << (this->usePeriodicBC ? "true" : "false") << std::endl;
    out << "calculateEigenvectors : " << (this->calculateEigenvectors ? "true" : "false") << std::endl;
    out << "saveEigenenergies     : " << (this->saveEigenenergies ? "true" : "false") << std::endl;
    out << "saveEigenstates       : " << (this->saveEigenstates ? "true" : "false") << std::endl;
    out << "storeFormat           : " << this->storeFormat << std::endl;
    out << "from                  : " << this->from << std::endl;
    out << "to                    : " << this->to << std::endl;
    out << "totalSimulations      : " << this->totalSimulations << std::endl;
    out << "seed                  : " << this->seed << std::endl;
    out << "splitWorkload         : " << (this->splitWorkload ? "true" : "false") << std::endl;
    out << "secureSimulationState : " << (this->secureSimulationState ? "true" : "false") << std::endl;
}

void Parameters::printHamiltonianTerms(std::ostream &out) const {
    out << "Hamiltonian terms:" << std::endl;
    for (const auto &term : this->hamiltonianTerms) {
        out << "-- " << term.first << std::endl;
        out << term.second;
    }
}

void Parameters::print(std::ostream &out) const {
    this->printGeneral(out);
    out << std::endl;
    this->printHamiltonianTerms(out);
}

std::string Parameters::getByName(const std::string &name) const {
    // General parameters
    if (name == "N")
        return this->doubleToString(this->N);
    else if (name == "K")
        return this->doubleToString(this->K);
    else if (name == "averagingModel")
        return this->averagingModel;
    else if (name == "usePeriodicBC")
        return this->usePeriodicBC ? "true" : "false";
    else if (name == "calculateEigenvectors")
        return this->calculateEigenvectors ? "true" : "false";
    else if (name == "saveEigenenergies")
        return this->saveEigenenergies ? "true" : "false";
    else if (name == "saveEigenstates")
        return this->saveEigenstates ? "true" : "false";
    else if (name == "from")
        return std::to_string(this->from);
    else if (name == "to")
        return std::to_string(this->to);
    else if (name == "totalSimulations")
        return std::to_string(this->totalSimulations);
    else if (name == "seed")
        return std::to_string(this->seed);
    else if (name == "splitWorkload")
        return std::to_string(this->splitWorkload);
    else if (name == "secureSimulationState")
        return std::to_string(this->secureSimulationState);

    // Hamiltonian term parameters
    for (auto &term : this->hamiltonianTerms) {
        const auto &termParams = term.second;
        if (termParams.hasField(name))
            return termParams.getString(name);
    }

    throw ParametersParseException("Unknown parameter " + name);
}

std::string Parameters::doubleToString(double d) const {
    std::ostringstream ostream;
    ostream << d;
    return ostream.str();
}

std::string Parameters::getOutputFileSignature() const {
    std::ostringstream filename;
    filename << "N." << this->N << "_K." << this->K;
    this->appendHamiltonianTermsSignature(filename);
    return filename.str();
}

std::string Parameters::getOutputFileSignatureWithRange() const {
    std::ostringstream filename;
    filename << "N." << this->N << "_K." << this->K << "_from." << this->from << "_to." << this->to;
    this->appendHamiltonianTermsSignature(filename);
    return filename.str();
}

void Parameters::appendHamiltonianTermsSignature(std::ostringstream &filename) const {
    for (const auto &term : hamiltonianTerms) {
        const auto &termParams = term.second;
        for (const auto &termParam : termParams.getKeys()) {
            std::string paramValue = termParams.getString(termParam);

            // Erase beginnings of paths ~/, ./ and ../ (possibly repeated)
            // ' ' and '/' should be dots
            // Eventually "pathParam = ../../subdir/file.dat" becomes nice - "[...]_pathParam.subdir.file.dat_[...]"
            paramValue = replaceAll(paramValue, "../", "");
            paramValue = replaceAll(paramValue, "./", "");
            paramValue = replaceAll(paramValue, "~/", "");
            paramValue = replaceAll(paramValue, " ", ".");
            paramValue = replaceAll(paramValue, "/", ".");

            filename << "_" << termParam << "." << paramValue;
        }
    }
}

bool Parameters::hasParam(const std::string &name) const {
    try {
        static_cast<void>(this->getByName(name));
    } catch (const ParametersParseException &e) {
        return false;
    }
    return true;
}
