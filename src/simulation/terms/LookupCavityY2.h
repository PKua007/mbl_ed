//
// Created by Piotr Kubala on 14/02/2020.
//

#ifndef MBL_ED_LOOKUPCAVITYY2_H
#define MBL_ED_LOOKUPCAVITYY2_H

#include <utility>

#include "simulation/CavityConstants.h"
#include "simulation/DoubleHoppingTerm.h"
#include "utils/Assertions.h"

/**
 * @brief Higher order term when computing cavity mediated interactions. It consists of two one-site hops.
 * @details
 * <p>It reads as
 * \f[ -\frac{U_1}{K} \hat{Y}^2, \f]
 * where
 * \f[ \hat{Y} = \sum_{i=1}^K y_i(\hat{b}^\dagger_i\hat{b}_{i+1} + \text{c. c.}), \f]
 * \f[ y_i = \int dx \cos(kx)w_i(x)w_{i+1}(x), \f]
 * \f$ i \f$ is the number of site, \f$ K \f$ is the total number of sites, \f$ \hat{b}_i \f$ is annihilation
 * operator, \f$ w_i \f$ is Wannier function localized in site \f$ i \f$ and \f$ U_1 \f$ is passed in the constructor.
 * <p> The constant are taken from the lookup table CavityConstants from the constructor.
 */
class LookupCavityY2 : public DoubleHoppingTerm {
private:
    double U1{};
    CavityConstants cavityConstants;
    CavityConstants::Realisation currentRealisation;

public:
    LookupCavityY2(double U1, CavityConstants cavityConstants, std::size_t realisationIndex = 0)
            : U1{U1}, cavityConstants{std::move(cavityConstants)}
    {
        Expects(U1 >= 0);
        Expects(realisationIndex < this->cavityConstants.size());
        this->currentRealisation = this->cavityConstants[realisationIndex];
    }

    /**
     * @brief Calculates the hop matrixElement for two hops. CavityConstants from the constructor must have enough sites
     * defined.
     */
    double calculate(const HopData &firstHopData, const HopData &secondHopData,
                     const HamiltonianGenerator &generator) override;

    /**
     * @brief Changes the realisation, i.e. phi0 and wanniers, to the one pointed by @a index in CavityConstants from
     * the constructor.
     */
    void changeRealisation(std::size_t index);
};


#endif //MBL_ED_LOOKUPCAVITYY2_H
