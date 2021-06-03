//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_HOPPINGTERM_H
#define MBL_ED_HOPPINGTERM_H

#include "FockBasis.h"

class HamiltonianGenerator;
struct HopData;

/**
 * @brief A class representing hopping term (off-diagonal term) of second-quantized hamiltonian.
 */
class HoppingTerm {
public:
    virtual ~HoppingTerm() = default;

    /**
     * @brief It is supposed to return a constant \f$ Y \f$ from a term
     * \f$ Y\ \hat{b}^\dagger_\text{hopData.toSite} \hat{b}_\text{hopData.fromSite} \f$ when acting on vector
     * \f$|\text{hopData.fromVector}\rangle \f$ and moving it to \f$ |\text{hopData.toVector}\rangle \f$.
     * @details This constant may depend on all of these things: from which site to which there is a hop and what are
     * initial and final Fock vectors (eg. density dependent hoppings). Note, that factors which appear just from ladder
     * operators should not be included, only this additional \f$ Y \f$ factor.
     */
    [[nodiscard]] virtual double calculate(const HopData &hopData, const HamiltonianGenerator &generator) const = 0;

    /**
     * @brief Returns a vector of all hopping distances with nonzero factors in this HoppingTerm that should be sampled.
     */
    [[nodiscard]] virtual std::vector<std::size_t> getHoppingDistances() const = 0;
};

#endif //MBL_ED_HOPPINGTERM_H
