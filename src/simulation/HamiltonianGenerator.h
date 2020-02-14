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

struct HopData {
    std::size_t fromSite{};
    std::size_t toSite{};
    FockBase::Vector fromVector{};
    FockBase::Vector toVector{};
    double ladderConstant{};
};

/**
 * @brief Hamiltonian generator, which can accept multiple DiagonalTerm -s and HoppingTerm -s.
 */
class HamiltonianGenerator {
private:
    const bool usePBC;
    std::unique_ptr<FockBase> fockBase;
    std::vector<std::unique_ptr<DiagonalTerm>> diagonalTerms;
    std::vector<std::unique_ptr<HoppingTerm>> hoppingTerms;

    [[nodiscard]] std::optional<HopData>
    hoppingAction(const FockBase::Vector &fromVector, std::size_t fromSite, std::size_t toSite) const;

public:
    HamiltonianGenerator(std::unique_ptr<FockBase> fockBase, bool usePBC)
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
    [[nodiscard]] size_t getSiteDistance(size_t fromSiteIndex, size_t toSiteIndex) const;

    void addDiagonalTerm(std::unique_ptr<DiagonalTerm> term);
    void addHoppingTerm(std::unique_ptr<HoppingTerm> term);

    /**
     * @brief Returns modyfiable list of diagonal terms.
     */
    [[nodiscard]] std::vector<std::unique_ptr<DiagonalTerm>> &getDiagonalTerms();

    /**
     * @brief Returns modyfiable list of diagonal terms.
     */
    [[nodiscard]] std::vector<std::unique_ptr<HoppingTerm>> &getHoppingTerms();

    [[nodiscard]] const FockBase &getFockBase() const;

    [[nodiscard]] bool usingPBC() const { return this->usePBC; }
};


#endif //MBL_ED_HAMILTONIANGENERATOR_H
