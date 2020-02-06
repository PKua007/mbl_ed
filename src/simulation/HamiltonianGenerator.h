//
// Created by pkua on 01.11.2019.
//

#ifndef MBL_ED_HAMILTONIANGENERATOR_H
#define MBL_ED_HAMILTONIANGENERATOR_H

#include <armadillo>
#include <random>
#include <memory>

#include "FockBase.h"

/**
 * @brief A base class for hamiltonians of a fixed number of particles in a linear optical lattice.
 * @details Is allows derived classes to provide their own diagonal terms for each fock state and kinetic constants for
 * each hopping term while doing the rest of the job of producing hamiltonian matrix.
 */
class HamiltonianGenerator {
private:
    [[nodiscard]] std::optional<std::pair<FockBase::Vector, double>>
    hoppingAction(const FockBase::Vector &vector, std::size_t fromSiteIndex, std::size_t toSiteIndex) const;

protected:
    /**
     * @brief If true, PBC are used, OBC otherwise.
     */
    const bool usePBC;
    std::unique_ptr<FockBase> fockBase;

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

    /**
     * @brief Derived classes should implement this method to provide own diagonal terms for all fock states.
     */
    [[nodiscard]] virtual double getDiagonalElement(const FockBase::Vector &vector) const = 0;

    /**
     * @brief Derived classes should implement this method to provide hopping constants for each possible hop.
     * @details Note, that hops of any lengths will be sampled and the method should return 0 if some of them are
     * forbidden.
     */
    [[nodiscard]] virtual double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const = 0;
};


#endif //MBL_ED_HAMILTONIANGENERATOR_H
