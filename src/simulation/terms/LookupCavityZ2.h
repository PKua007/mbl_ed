//
// Created by Piotr Kubala on 10/02/2020.
//

#ifndef MBL_ED_LOOKUPCAVITYZ2_H
#define MBL_ED_LOOKUPCAVITYZ2_H


#include <utility>

#include "simulation/CavityConstants.h"
#include "simulation/DiagonalTerm.h"
#include "utils/Assertions.h"

/**
 * @brief Like CavityLongInteraction, but with lookup values of cosine (or rather Wannier integral in this case)
 * @details
 * <p> Wannier integral: \f$ \int dx \cos(kx) w_i(x)^2 \f$ which, for perfectly localized Wanniers \f$ w_i \f$
 * takes the form from CavityLongInteraction \f$ \cos(2\pi\beta i + \phi_0) \f$.
 * <p> The values are taken from CavityConstants table.
 */
class LookupCavityZ2 : public DiagonalTerm {
private:
    double U1{};
    CavityConstants cavityConstants;
    CavityConstants::Realisation currentRealisation;

public:
    LookupCavityZ2(double U1, CavityConstants cavityConstants, std::size_t realisationIndex = 0)
            : U1{U1}, cavityConstants{std::move(cavityConstants)}
    {
        Expects(U1 >= 0);
        Expects(realisationIndex < this->cavityConstants.size());
        this->currentRealisation = this->cavityConstants[realisationIndex];
    }

    /**
     * @brief Calculates the diagonal elements for @a vector. CavityConstants from the constructor must have enough
     * sites defined.
     */
    double calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) override;

    /**
     * @brief Changes the realisation, i.e. phi0 and wanniers, to the one pointed by @a index in CavityConstants from
     * the constructor.
     */
    void changeRealisation(std::size_t index);
};


#endif //MBL_ED_LOOKUPCAVITYZ2_H
