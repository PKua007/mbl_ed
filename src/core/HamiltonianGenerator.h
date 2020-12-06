//
// Created by pkua on 01.11.2019.
//

#ifndef MBL_ED_HAMILTONIANGENERATOR_H
#define MBL_ED_HAMILTONIANGENERATOR_H

#include <armadillo>
#include <random>
#include <memory>

#include "FockBasis.h"
#include "DiagonalTerm.h"
#include "HoppingTerm.h"
#include "DoubleHoppingTerm.h"
#include "Eigensystem.h"

/**
 * @brief Struct representing a hop between two sites.
 */
struct HopData {
    std::size_t fromSite{};
    std::size_t toSite{};
    FockBasis::Vector fromVector{};
    FockBasis::Vector toVector{};

    /**
     * @brief The constant given by acting with \f$ \hat{b}_\text{toSite} \hat{b}_\text{fromSite} \f$ on
     * \f$ |\text{fromVector}> \f$.
     */
    double ladderConstant{};
};

/**
 * @brief Hamiltonian generator, which can accept multiple DiagonalTerm -s, HoppingTerm -s and DoubleHoppingTerm -s.
 */
class HamiltonianGenerator {
private:
    const bool usePBC;
    std::shared_ptr<const FockBasis> fockBasis;
    std::vector<std::shared_ptr<DiagonalTerm>> diagonalTerms;
    std::vector<std::shared_ptr<HoppingTerm>> hoppingTerms;
    std::vector<std::shared_ptr<DoubleHoppingTerm>> doubleHoppingTerms;

    [[nodiscard]] std::optional<HopData>
    hoppingAction(const FockBasis::Vector &fromVector, std::size_t fromSite, std::size_t toSite) const;
    [[nodiscard]] auto calculateDoubleHopMatrixElement(const HopData &firstHop, const HopData &secondHop) const;
    void performSecondHop(arma::sp_mat &result, std::size_t fromIdx, const HopData &firstHop) const;
    void addDiagonalTerms(arma::sp_mat &result, std::size_t vectorIdx) const;
    void addHoppingTerms(arma::sp_mat &result, std::size_t fromIdx) const;
    void addDoubleHoppingTerms(arma::sp_mat &result, std::size_t fromIdx) const;

public:
    HamiltonianGenerator(std::shared_ptr<const FockBasis> fockBasis, bool usePBC)
            : usePBC{usePBC}, fockBasis{std::move(fockBasis)}
    { }
    virtual ~HamiltonianGenerator() = default;

    /**
     * @brief Generates the hamiltonian matrix - it uses all provided DiagonalTerm -s and HoppingTerms -s.
     * @details For HoppingTerms -s, it produces them by doing one-site hop on each basis vector and finding which
     * another basis vector is obtained that way.
     */
    [[nodiscard]] arma::sp_mat generate() const;

    /**
     * @brief Generates hamiltonian and diagonalizes it. It is not dumb, if hamiltonian is diagonal it doesn't
     * invoke diagonalization routines.
     * @return Eigensystem with or without eigenvectors, depending on @a calculateEigenvectors flag.
     */
    [[nodiscard]] Eigensystem calculateEigensystem(bool calculateEigenvectors) const;

    /**
     * @brief Returns the distance between sites of given indices.
     * @details Note, that when used with PBC, the shorter distance will be returned.
     */
    [[nodiscard]] size_t getSiteDistance(std::size_t fromSite, std::size_t toIdx) const;

    void addDiagonalTerm(std::shared_ptr<DiagonalTerm> term);
    void addHoppingTerm(std::shared_ptr<HoppingTerm> term);
    void addDoubleHoppingTerm(std::shared_ptr<DoubleHoppingTerm> term);

    /**
     * @brief Returns modyfiable list of diagonal terms.
     */
    [[nodiscard]] const std::vector<std::shared_ptr<DiagonalTerm>> &getDiagonalTerms();

    /**
     * @brief Returns modyfiable list of hopping terms.
     */
    [[nodiscard]] const std::vector<std::shared_ptr<HoppingTerm>> &getHoppingTerms();

    /**
     * @brief Returns modyfiable list of double hopping terms.
     */
    [[nodiscard]] const std::vector<std::shared_ptr<DoubleHoppingTerm>> &getDoubleHoppingTerms();

    [[nodiscard]] const std::shared_ptr<const FockBasis> &getFockBasis() const { return this->fockBasis; };
    [[nodiscard]] bool usingPBC() const { return this->usePBC; }
};


#endif //MBL_ED_HAMILTONIANGENERATOR_H
