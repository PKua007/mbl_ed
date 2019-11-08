//
// Created by pkua on 05.11.2019.
//

#ifndef MBL_ED_CAVITYHAMILTONIANGENERATOR_H
#define MBL_ED_CAVITYHAMILTONIANGENERATOR_H

#include <random>

#include "Assertions.h"
#include "HamiltonianGenerator.h"

template<typename DisorderGenerator>
class CavityHamiltonianGenerator : public HamiltonianGenerator {
private:
    double J{};
    double U{};
    double U1{};
    std::vector<double> onsiteEnergies;

public:
    CavityHamiltonianGenerator(const FockBase &fockBase, double J, double U, double U1,
                               DisorderGenerator &&disorderGenerator, bool usePbc = true)
            : HamiltonianGenerator(fockBase, usePbc), J(J), U(U), U1(U1)
    {
        this->onsiteEnergies.resize(fockBase.getNumberOfSites());
        std::generate(this->onsiteEnergies.begin(), this->onsiteEnergies.end(), disorderGenerator);
    }

    [[nodiscard]] double getDiagonalElement(const FockBase::Vector &vector) const override {
        std::vector<double> elementwiseEnergies;
        elementwiseEnergies.reserve(vector.size());
        std::transform(vector.begin(), vector.end(), this->onsiteEnergies.begin(), elementwiseEnergies.begin(),
                       std::multiplies<>());
        double onsiteEnergy = std::accumulate(elementwiseEnergies.begin(), elementwiseEnergies.end(), 1.,
                                              std::multiplies<>());

        auto bosonInteractionAccumulator = [](auto sum, auto numberOfParticles) {
            return sum + numberOfParticles*(numberOfParticles - 1);
        };
        double shortInteractionEnergy = U * std::accumulate(vector.begin(), vector.end(), 0., bosonInteractionAccumulator);

        int multiplier = -1;
        auto plusMinusAccumulator = [&multiplier](auto sum, auto element) {
            return sum + (multiplier *= -1) * element;
        };
        double populationImbalance = std::accumulate(vector.begin(), vector.end(), 0., plusMinusAccumulator);
        double longInteractionEnergy = -this->U1/vector.size() * populationImbalance*populationImbalance;

        return onsiteEnergy + shortInteractionEnergy + longInteractionEnergy;
    }

    [[nodiscard]] double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const override {
        Expects(this->getSiteDistance(fromSiteIndex, toSiteIndex) == 1);

        return -this->J;
    }
};

#endif //MBL_ED_CAVITYHAMILTONIANGENERATOR_H
