//
// Created by pkua on 01.11.2019.
//

#ifndef MBL_ED_HAMILTONIANGENERATOR_H
#define MBL_ED_HAMILTONIANGENERATOR_H

#include <armadillo>
#include <random>
#include <memory>

#include "FockBase.h"

class HamiltonianGenerator {
private:
    [[nodiscard]] std::optional<std::pair<FockBase::Vector, double>>
    hoppingAction(const FockBase::Vector &vector, std::size_t fromSiteIndex, std::size_t toSiteIndex) const;

protected:
    const bool usePBC;
    std::unique_ptr<FockBase> fockBase;

public:
    explicit HamiltonianGenerator(std::unique_ptr<FockBase> fockBase, bool usePBC)
            : usePBC{usePBC}, fockBase{std::move(fockBase)}
    { }
    virtual ~HamiltonianGenerator() = default;

    [[nodiscard]] arma::mat generate() const;
    [[nodiscard]] size_t getSiteDistance(size_t fromSiteIndex, size_t toSiteIndex) const;

    [[nodiscard]] virtual double getDiagonalElement(const FockBase::Vector &vector) const = 0;
    [[nodiscard]] virtual double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const = 0;
};


#endif //MBL_ED_HAMILTONIANGENERATOR_H
