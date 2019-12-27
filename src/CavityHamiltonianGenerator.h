//
// Created by pkua on 05.11.2019.
//

#ifndef MBL_ED_CAVITYHAMILTONIANGENERATOR_H
#define MBL_ED_CAVITYHAMILTONIANGENERATOR_H

#include <random>
#include <memory>
#include <iterator>

#include "Assertions.h"
#include "HamiltonianGenerator.h"

/**
 * @brief Generator of boson hamiltonian from https://arxiv.org/pdf/1902.00357.pdf
 * @tparam DisorderGenerator class whose operator() returns random energies as double
 */
template<typename DisorderGenerator>
class CavityHamiltonianGenerator : public HamiltonianGenerator {
private:
    double J{};
    double U{};
    double U1{};
    std::vector<double> onsiteEnergies;
    std::unique_ptr<DisorderGenerator> disorderGenerator;

    /* Sum of E_j n_j, where j = 0, ..., [number of sites]; n_j - number of particles in j-th site and E_j - random site
     * energies from the constructor */
    [[nodiscard]] double getOnsiteEnergy(const FockBase::Vector &vector) const {
        std::vector<double> elementwiseEnergies;
        elementwiseEnergies.reserve(vector.size());
        std::transform(vector.begin(), vector.end(), this->onsiteEnergies.begin(),
                       std::back_inserter(elementwiseEnergies), std::multiplies<>());
        return std::accumulate(elementwiseEnergies.begin(), elementwiseEnergies.end(), 0., std::plus<>());
    }

    /* Sum of U * n_j * (n_j - 1), where j = 1, ..., [number of sites]; n_j - number of particles in j-th site */
    [[nodiscard]] double getShortInteractionEnergy(const FockBase::Vector &vector) const {
        auto bosonAccumulator = [](auto sum, auto numberOfParticles) {
            return sum + numberOfParticles*(numberOfParticles - 1);
        };
        return this->U * std::accumulate(vector.begin(), vector.end(), 0., bosonAccumulator);
    }

    /* -U1/[number of sites] * (sum of (-1)^j n_j)^2, where j = 1, ..., [number of sites]; n_j - number of particles in
     * j-th site */
    [[nodiscard]] double getLongInteractionEnergy(const FockBase::Vector &vector) const {
        int multiplier = -1;
        auto plusMinusAccumulator = [&multiplier](auto sum, auto element) {
            return sum + (multiplier *= -1) * element;
        };
        double populationImbalance = std::accumulate(vector.begin(), vector.end(), 0., plusMinusAccumulator);
        return -this->U1 / vector.size() * populationImbalance * populationImbalance;
    }

public:
    CavityHamiltonianGenerator(const FockBase &fockBase, double J, double U, double U1,
                               std::unique_ptr<DisorderGenerator> disorderGenerator, bool usePbc = true)
            : HamiltonianGenerator(fockBase, usePbc), J(J), U(U), U1(U1),
              disorderGenerator(std::move(disorderGenerator))
    {
        this->resampleOnsiteEnergies();
    }

    void resampleOnsiteEnergies() {
        this->onsiteEnergies.resize(this->fockBase.getNumberOfSites());
        std::generate(this->onsiteEnergies.begin(), this->onsiteEnergies.end(), std::ref(*(this->disorderGenerator)));
    }

    [[nodiscard]] const std::vector<double> &getOnsiteEnergies() const {
        return this->onsiteEnergies;
    }

    [[nodiscard]] double getDiagonalElement(const FockBase::Vector &vector) const override {
        return getOnsiteEnergy(vector) + getShortInteractionEnergy(vector) + getLongInteractionEnergy(vector);
    }

    [[nodiscard]] double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const override {
        Expects(this->getSiteDistance(fromSiteIndex, toSiteIndex) == 1);

        return -this->J;
    }

    [[nodiscard]] std::string fileSignature() const {
        std::ostringstream filename;
        filename << "J." << this->J << "_U." << this->U << "_U1." << this->U1;
        filename << "_N." << this->fockBase.getNumberOfParticles() << "_K." << this->fockBase.getNumberOfSites();
        return filename.str();
    }
};

#endif //MBL_ED_CAVITYHAMILTONIANGENERATOR_H
