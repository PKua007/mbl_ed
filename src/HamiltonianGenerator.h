//
// Created by pkua on 01.11.2019.
//

#ifndef MBL_ED_HAMILTONIANGENERATOR_H
#define MBL_ED_HAMILTONIANGENERATOR_H

#include <armadillo>

#include "FockBase.h"

class HamiltonianGenerator {
private:
    const bool usePBC;

    [[nodiscard]] std::optional<std::pair<FockBase::Vector, double>>
    hoppingAction(const FockBase::Vector &vector, std::size_t fromSiteIndex, std::size_t toSiteIndex) const;

protected:
    const FockBase &fockBase;

public:
    explicit HamiltonianGenerator(const FockBase &fockBase, bool usePBC = true) : usePBC{usePBC}, fockBase{fockBase} { }

    [[nodiscard]] arma::mat generate() const;

    [[nodiscard]] virtual double getDiagonalElement(const FockBase::Vector &vector) const = 0;
    [[nodiscard]] virtual double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const = 0;
};


#endif //MBL_ED_HAMILTONIANGENERATOR_H
