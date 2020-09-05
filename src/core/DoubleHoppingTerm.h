//
// Created by Piotr Kubala on 14/02/2020.
//

#ifndef MBL_ED_DOUBLEHOPPINGTERM_H
#define MBL_ED_DOUBLEHOPPINGTERM_H

class HamiltonianGenerator;
class HopData;

/**
 * @brief A class representing double hopping term (two \f$ b^\dagger b \f$ pairs) of second-quantized hamiltonian.
 */
class DoubleHoppingTerm {
public:
    virtual ~DoubleHoppingTerm() = default;

   /**
     * @brief It is supposed to return a constant \f$ Y \f$ from a term
     * \f$ Y\ \hat{b}^\dagger_\text{secondHopData.toSite} \hat{b}_\text{secondHopData.fromSite}
     *        \hat{b}^\dagger_\text{firstHopData.toSite} \hat{b}_\text{firstHopData.fromSite}\f$ when acting on vector
     * \f$|\text{firstHopData.fromVector}\rangle \f$ and moving it to \f$ |\text{secondHopData.toVector}\rangle \f$.
     * @details This constant may depend on all sites' occupance during each step of hopping. Note, that factors which
     * appear just from ladder operators should not be included, only this additional \f$ Y \f$ factor.
     */
    virtual double calculate(const HopData &firstHopData, const HopData &secondHopData,
                             const HamiltonianGenerator &generator) = 0;
};

#endif //MBL_ED_DOUBLEHOPPINGTERM_H
