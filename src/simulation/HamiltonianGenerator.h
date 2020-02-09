//
// Created by pkua on 01.11.2019.
//

#ifndef MBL_ED_HAMILTONIANGENERATOR_H
#define MBL_ED_HAMILTONIANGENERATOR_H

#include <armadillo>
#include <random>
#include <memory>

#include "FockBase.h"
#include "DiagonalTerm.h"
#include "HoppingTerm.h"

/**
 * @brief A base class for hamiltonians of a fixed number of particles in a linear optical lattice.
 * @details Is allows derived classes to provide their own diagonal terms for each fock state and kinetic constants for
 * each hopping term while doing the rest of the job of producing hamiltonian matrix.
 */
class HamiltonianGenerator {
private:
    const bool usePBC;
    std::unique_ptr<FockBase> fockBase;
    std::vector<std::unique_ptr<DiagonalTerm>> diagonalTerms;
    std::vector<std::unique_ptr<HoppingTerm>> hoppingTerms;

    [[nodiscard]] std::optional<std::pair<FockBase::Vector, double>>
    hoppingAction(const FockBase::Vector &vector, std::size_t fromSiteIndex, std::size_t toSiteIndex) const;

public:
    HamiltonianGenerator(std::unique_ptr<FockBase> fockBase, bool usePBC)
            : usePBC{usePBC}, fockBase{std::move(fockBase)}
    { }
    virtual ~HamiltonianGenerator() = default;

    /**
     * @brief Using derived class' provided getDiagonalElement() and getHoppingTerm() methods, generates the complete
     * hamiltonian matrix and returns it.
     */
    [[nodiscard]] arma::mat generate() const;

    /**
     * @brief Returns the distance between sites of given indices.
     * @details Note, that when used with PBC, the shorter distance will be returned.
     */
    [[nodiscard]] size_t getSiteDistance(size_t fromSiteIndex, size_t toSiteIndex) const;

    void addDiagonalTerm(std::unique_ptr<DiagonalTerm> term);
    void addHoppingTerm(std::unique_ptr<HoppingTerm> term);
    [[nodiscard]] const std::vector<std::unique_ptr<DiagonalTerm>> &getDiagonalTerms() const;
    [[nodiscard]] const std::vector<std::unique_ptr<HoppingTerm>> &getHoppingTerms() const;
    [[nodiscard]] const FockBase &getFockBase() const;
};


#endif //MBL_ED_HAMILTONIANGENERATOR_H
