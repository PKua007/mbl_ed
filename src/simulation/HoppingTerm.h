//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_HOPPINGTERM_H
#define MBL_ED_HOPPINGTERM_H

#include "FockBase.h"

class HamiltonianGenerator;

/**
 * @brief A class representing hopping term (off-diagonal term) of second-quantized hamiltonian.
 */
class HoppingTerm {
public:
    virtual ~HoppingTerm() = default;

    /**
     * @brief It is supposed to return a constant \f$ Y \f$ from a term
     * \f$ Y\ \hat{b}^\dagger_\text{toSite} \hat{b}_\text{fromSite} \f$ when acting on vector \f$|\text{from}\rangle \f$
     * and moving it to \f$ |\text{to}\rangle \f$.
     * @details This constant may depend on all of there things: from which site to which there is a hop and what are\
     * initial and final Fock vectors. Note, that factors which appear just from ladder operators should not be
     * included, only this additional \f$ Y \f$ factor.
     */
    virtual double calculate(const FockBase::Vector &from, const FockBase::Vector &to, std::size_t fromSite,
                             std::size_t toSite, const HamiltonianGenerator &generator) = 0;
};

#endif //MBL_ED_HOPPINGTERM_H
