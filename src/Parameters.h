//
// Created by Piotr Kubala on 28/12/2019.
//

#ifndef MBL_ED_PARAMETERS_H
#define MBL_ED_PARAMETERS_H

#include <iosfwd>

class Parameters {
private:
    void validate() const;

public:
    std::size_t numberOfSites{};
    std::size_t numberOfBosons{};
    double J{};
    double W{};
    double U{};
    double U1{};
    double beta{};
    double phi0{};
    bool usePeriodicBC{};
    std::size_t numberOfSimulations{};
    std::size_t seed{};

    explicit Parameters(std::istream &input);

    void print(std::ostream &out) const;
};


#endif //MBL_ED_PARAMETERS_H
