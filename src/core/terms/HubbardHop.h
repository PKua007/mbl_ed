//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_HUBBARDHOP_H
#define MBL_ED_HUBBARDHOP_H


#include "utils/Assertions.h"
#include "core/HoppingTerm.h"

/**
 * @brief A standard one-site hop from Hubbard Hamiltonian.
 * @details It is defined as
 * \f[ -J_j \sum_{i=1}^K ( \hat{b}_i^\dagger \hat{b}_{i+j} + \text{c. c.} ), \f]
 *
 * where \f$ i \f$ is the number of site, \f$ K \f$ is the total number of sites, \f$ \hat{b}_i \f$ is annihilation
 * operator and the constant \f$ J \f$ is passed in the constructor. Moreover \f$ j \f$ may take one or several values
 * representing hoppings by that many sites.
 */
class HubbardHop : public HoppingTerm {
private:
    std::vector<double> Js{};
    std::vector<std::size_t> hoppingDistances{};

public:
    /**
     * @brief Hoppings by 1 site of amplitude @a J
     */
    explicit HubbardHop(double J);

    /**
     * @brief Hoppings by 1, 2, 3, ..., n sites of amplitudes @a Js[1], ... @a Js[n]
     */
    explicit HubbardHop(std::vector<double> Js);

    /**
     * @brief Hoppings by explicitly states numbers of sites with corresponding amplitudes
     */
    HubbardHop(std::vector<std::size_t> hoppingDistances, std::vector<double> Js);

    /**
     * @brief Returns constant -J for a hop between neighbouring sites (according to PBC or OBC).
     */
    double calculate(const HopData &hopData, const HamiltonianGenerator &generator) const override;

    [[nodiscard]] std::vector<std::size_t> getHoppingDistances() const override { return this->hoppingDistances; }
};

#endif //MBL_ED_HUBBARDHOP_H
