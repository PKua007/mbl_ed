//
// Created by Piotr Kubala on 28/12/2019.
//

#ifndef MBL_ED_PARAMETERS_H
#define MBL_ED_PARAMETERS_H

#include <iosfwd>
#include <exception>
#include <map>

#include "utils/Config.h"

/**
 * @brief Exception thrown when parsing of accessing unknown parameter.
 */
struct ParametersParseException : public std::runtime_error {
public:
    explicit ParametersParseException(const std::string &what) : std::runtime_error(what) { }
};

/**
 * @brief A class which parses and stores parameters of the simulation.
 * @details The description of parameters can be found in the sample input.ini in the root folder
 */
class Parameters {
private:
    void autocompleteAndValidate();
    [[nodiscard]] std::string doubleToString(double d) const;
    void appendHamiltonianTermsSignature(std::ostringstream &filename) const;

public:
    /* All of these are described in input.ini */

    std::size_t N{};
    std::size_t K{};
    bool usePeriodicBC{};
    bool calculateEigenvectors = true;
    bool saveEigenenergies = false;
    std::string averagingModel{};
    std::size_t to{};
    std::size_t from = 0;
    std::size_t totalSimulations = 1;
    std::size_t seed{};

    /**
     * @brief All keys from sections @a [term.termName] are mapped to separate config under @a termName key in the map.
     */
    std::map<std::string, Config> hamiltonianTerms;

    Parameters() = default;

    /**
     * @brief Parses Config type input from @a input stream.
     */
    explicit Parameters(std::istream &input);

    /**
     * @brief Prints the summary of general parameters to the @a out stream.
     */
    void printGeneral(std::ostream &out) const;

    /**
     * @brief Prints the summary of hamiltonian terms to the @a out stream.
     */
    void printHamiltonianTerms(std::ostream &out) const;

    /**
     * @brief Prints the summary of both general parameters and hamiltonian terms to the @a out stream.
     */
    void print(std::ostream &out) const;

    /**
     * @brief Returns general or term parameters by name. Name should not contain the section.
     */
    [[nodiscard]] std::string getByName(const std::string &name) const;

    /**
     * @brief Check is the paramm of given name exists, Name should not contain the section.
     */
    [[nodiscard]] bool hasParam(const std::string &name) const;

    /**
     * @brief Prepares prefix for output files which gives some info about parameters.
     * @details It prints: N, K and all parameters of hamiltonian terms, without the name of the sections. The format
     * may change often, consult the source code.
     */
    [[nodiscard]] std::string getOutputFileSignature() const;

    /**
     * @brief The same as getOutputFileSignature(), but parameters Parameters::from and Parameters::to are included.
     */
    [[nodiscard]] std::string getOutputFileSignatureWithRange() const;
};


#endif //MBL_ED_PARAMETERS_H
