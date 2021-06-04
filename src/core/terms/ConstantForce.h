//
// Created by pkua on 04.06.2021.
//

#ifndef MBL_ED_CONSTANTFORCE_H
#define MBL_ED_CONSTANTFORCE_H

#include "core/DiagonalTerm.h"

/**
 * @brief A class representing constant force (linear potential).
 * @details It is realised by linearly increasing onsite potential. It is shifted in a way, that middle site has energy
 * zero and two ends of a chain have opposite energies.
 */
class ConstantForce : public DiagonalTerm {
private:
    double F{};

public:
    /**
     * @brief Constructs a class for (force) * (site distance) = @a F, namely potential difference between neighbouring
     * sites will be @a F.
     */
    explicit ConstantForce(double F);

    double calculate(const FockBasis::Vector &vector, const HamiltonianGenerator &generator) override;
};


#endif //MBL_ED_CONSTANTFORCE_H
