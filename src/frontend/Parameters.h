//
// Created by Piotr Kubala on 28/12/2019.
//

#ifndef MBL_ED_PARAMETERS_H
#define MBL_ED_PARAMETERS_H

#include <iosfwd>
#include <exception>

struct UnknownParameterException : public std::runtime_error {
public:
    explicit UnknownParameterException(const std::string &what) : std::runtime_error(what) { }
};

/**
 * @brief A class which parses and stores parameters of the simulation.
 * @details The description of parameters can be found in the sample input.txt in the root folder
 */
class Parameters {
private:
    void autocompleteAndValidate();
    [[nodiscard]] std::string doubleToString(double d) const;

public:
    /* All of these are described in input.txt */

    std::size_t numberOfSites{};
    std::size_t numberOfBosons{};
    double J{};
    double W{};
    double U{};
    double U1{};
    double beta{};
    std::string phi0;
    bool usePeriodicBC{};
    bool calculateEigenvectors = true;
    bool saveEigenenergies = false;
    std::size_t to{};
    std::size_t from = 0;
    std::size_t totalSimulations = 1;
    std::size_t seed{};

    Parameters() = default;

    /**
     * @brief Parses Config type input from @a input stream.
     */
    explicit Parameters(std::istream &input);

    /**
     * @brief Prints the summary of parameters to the @a out stream.
     */
    void print(std::ostream &out) const;

    [[nodiscard]] std::string getByName(const std::string &name) const;

    /**
     * @brief Prepares prefix for output files which gives some info about parameters.
     * @details See the source for the format, it may change often.
     */
    [[nodiscard]] std::string getOutputFileSignature() const;
};


#endif //MBL_ED_PARAMETERS_H
