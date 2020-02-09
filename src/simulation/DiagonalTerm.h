//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_DIAGONALTERM_H
#define MBL_ED_DIAGONALTERM_H

#include "FockBase.h"

class HamiltonianGenerator;

/**
 * @brief A class representing one diagonal term of second-quantized hamiltonian.
 */
class DiagonalTerm {
public:
    virtual ~DiagonalTerm() = default;

    /**
     * @brief Given fock base vector it is supposed to calculate its diagonal entry.
     * @details @a generator is passed in case it is needed.
     */
    virtual double calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) = 0;
};

#endif //MBL_ED_DIAGONALTERM_H
