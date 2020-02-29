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
#include "DoubleHoppingTerm.h"

/**
 * @brief Struct representing a hop between two sites.
 */
struct HopData {
    std::size_t fromSite{};
    std::size_t toSite{};
    FockBase::Vector fromVector{};
    FockBase::Vector toVector{};

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
    std::shared_ptr<const FockBase> fockBase;
    std::vector<std::unique_ptr<DiagonalTerm>> diagonalTerms;
    std::vector<std::unique_ptr<HoppingTerm>> hoppingTerms;
    std::vector<std::unique_ptr<DoubleHoppingTerm>> doubleHoppingTerms;

    [[nodiscard]] std::optional<HopData>
    hoppingAction(const FockBase::Vector &fromVector, std::size_t fromSite, std::size_t toSite) const;
    [[nodiscard]] auto calculateDoubleHopMatrixElement(const HopData &firstHop, const HopData &secondHop) const;
    void performSecondHop(arma::mat &result, std::size_t fromIdx, const HopData &firstHop) const;
    void addDiagonalTerms(arma::mat &result, std::size_t vectorIdx) const;
    void addHoppingTerms(arma::mat &result, std::size_t fromIdx) const;
    void addDoubleHoppingTerms(arma::mat &result, std::size_t fromIdx) const;

public:
    HamiltonianGenerator(std::shared_ptr<const FockBase> fockBase, bool usePBC)
            : usePBC{usePBC}, fockBase{std::move(fockBase)}
    { }
    virtual ~HamiltonianGenerator() = default;

    /**
     * @brief Generates the hamiltonian matrix - it uses all provided DiagonalTerm -s and HoppingTerms -s.
     * @details For HoppingTerms -s, it produces them by doing one-site hop on each base vector and finding which
     * another base vector is obtained that way.
     */
    [[nodiscard]] arma::mat generate() const;

    /**
     * @brief Returns the distance between sites of given indices.
     * @details Note, that when used with PBC, the shorter distance will be returned.
     */
    [[nodiscard]] size_t getSiteDistance(std::size_t fromSite, std::size_t toIdx) const;

    void addDiagonalTerm(std::unique_ptr<DiagonalTerm> term);
    void addHoppingTerm(std::unique_ptr<HoppingTerm> term);
    void addDoubleHoppingTerm(std::unique_ptr<DoubleHoppingTerm> term);

    /**
     * @brief Returns modyfiable list of diagonal terms.
     */
    [[nodiscard]] std::vector<std::unique_ptr<DiagonalTerm>> &getDiagonalTerms();

    /**
     * @brief Returns modyfiable list of hopping terms.
     */
    [[nodiscard]] std::vector<std::unique_ptr<HoppingTerm>> &getHoppingTerms();

    /**
     * @brief Returns modyfiable list of double hopping terms.
     */
    [[nodiscard]] std::vector<std::unique_ptr<DoubleHoppingTerm>> &getDoubleHoppingTerms();

    [[nodiscard]] const std::shared_ptr<const FockBase> &getFockBase() const { return this->fockBase; };
    [[nodiscard]] bool usingPBC() const { return this->usePBC; }
};


#endif //MBL_ED_HAMILTONIANGENERATOR_H
