//
// Created by pkua on 04.06.2021.
//

#include "ConstantForce.h"

#include "utils/Assertions.h"
#include "core/HamiltonianGenerator.h"

ConstantForce::ConstantForce(double F) : F{F} {
    Expects(F != 0);
}

double ConstantForce::calculate(const FockBasis::Vector &vector, const HamiltonianGenerator &generator) {
    Expects(!generator.usingPBC());

    double energyShift = -(static_cast<double>(vector.size()) - 1) / 2;
    int elementIndex{};
    auto energyAccumulator = [&elementIndex, energyShift](auto energy, auto occupation) {
        return energy + occupation * (energyShift + elementIndex++);
    };
    double constantForceEnergy = std::accumulate(vector.begin(), vector.end(), 0., energyAccumulator);
    return this->F * constantForceEnergy;
}
