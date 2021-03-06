//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_ONSITEDISORDER_H
#define MBL_ED_ONSITEDISORDER_H


#include <memory>
#include <algorithm>

#include "utils/Assertions.h"
#include "core/DiagonalTerm.h"
#include "core/RND.h"
#include "core/DisorderGenerator.h"

/**
 * @brief The class representing random potential in each site.
 * @details It is defined as
 * \f[ \sum_{i=1}^K E_i \hat{n}_i \f]
 *
 * where \f$ i \f$ is the number of site, \f$ K \f$ is the total number of sites, \f$ \hat{b}_i \f$ is annihilation
 * operator, \f$ \hat{n}_i = \hat{b}_i^\dagger\hat{b}_i \f$, \f$ E_i \f$ are random energies sampled by
 * @a DisorderGenerator. The disorder is sampled during construction and can be later resampled using
 * resampleOnsiteEnergies().
 */
class OnsiteDisorder : public DiagonalTerm {
private:
    std::vector<double> onsiteEnergies;
    std::unique_ptr<DisorderGenerator> disorderGenerator;

public:
    /**
     * @brief Contructs the class. @a numberOfSites is specified a priori and cannot be changed later. @a rnd
     * is used to do the initial sampling.
     */
    explicit OnsiteDisorder(std::unique_ptr<DisorderGenerator> disorderGenerator, std::size_t numberOfSites, RND &rnd);

    /**
     * @brief Clears old onsite disorder values and samples new ones. Used for each new simulation.
     */
    void resampleOnsiteEnergies(RND &rnd);

    double calculate(const FockBasis::Vector &vector, const HamiltonianGenerator &generator) override;
};


#endif //MBL_ED_ONSITEDISORDER_H
