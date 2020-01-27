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

class Parameters {
private:
    void autocompleteAndValidate();
    [[nodiscard]] std::string doubleToString(double d) const;

public:
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
    explicit Parameters(std::istream &input);

    void print(std::ostream &out) const;
    [[nodiscard]] std::string getByName(const std::string &name) const;
    [[nodiscard]] std::string getOutputFileSignature() const;
};


#endif //MBL_ED_PARAMETERS_H
