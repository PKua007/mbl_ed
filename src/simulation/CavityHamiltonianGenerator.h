//
// Created by pkua on 05.11.2019.
//

#ifndef MBL_ED_CAVITYHAMILTONIANGENERATOR_H
#define MBL_ED_CAVITYHAMILTONIANGENERATOR_H

#include <random>
#include <memory>
#include <iterator>

#include "utils/Assertions.h"
#include "HamiltonianGenerator.h"

struct CavityHamiltonianGeneratorParameters {
    double J{};
    double U{};
    double U1{};
    double beta{0.5};
    double phi0{};
};

/**
 * @brief Generator of boson hamiltonian from https://arxiv.org/pdf/1902.00357.pdf
 * @tparam DisorderGenerator class whose operator() returns random energies as double
 */
template<typename DisorderGenerator>
class CavityHamiltonianGenerator : public HamiltonianGenerator {
private:
    std::vector<double> onsiteEnergies;
    std::unique_ptr<DisorderGenerator> disorderGenerator;
    CavityHamiltonianGeneratorParameters params;

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
        return this->params.U / 2 * std::accumulate(vector.begin(), vector.end(), 0., bosonAccumulator);
    }

    /* -U1/[number of sites] * (sum of (-1)^j n_j)^2, where j = 1, ..., [number of sites]; n_j - number of particles in
     * j-th site */
    [[nodiscard]] double getLongInteractionEnergy(const FockBase::Vector &vector) const {
        std::size_t elementIndex{};
        auto plusMinusAccumulator = [&elementIndex, this](auto sum, auto element) {
            return sum + std::cos(2*M_PI*params.beta*(elementIndex++) + params.phi0) * element;
        };
        double populationImbalance = std::accumulate(vector.begin(), vector.end(), 0., plusMinusAccumulator);
        return -this->params.U1 / vector.size() * populationImbalance * populationImbalance;
    }

public:
    CavityHamiltonianGenerator(std::unique_ptr<FockBase> fockBase, CavityHamiltonianGeneratorParameters params,
                               std::unique_ptr<DisorderGenerator> disorderGenerator, bool usePbc = true)
            : HamiltonianGenerator(std::move(fockBase), usePbc), disorderGenerator(std::move(disorderGenerator)),
              params{params}
    {
        this->resampleOnsiteEnergies();
    }

    void resampleOnsiteEnergies() {
        this->onsiteEnergies.resize(this->fockBase->getNumberOfSites());
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

        return -this->params.J;
    }

    [[nodiscard]] std::string fileSignature() const {
        std::ostringstream filename;
        filename << "J." << this->params.J << "_U." << this->params.U << "_U1." << this->params.U1;
        filename << "_N." << this->fockBase->getNumberOfParticles() << "_K." << this->fockBase->getNumberOfSites();
        filename << "_beta." << this->params.beta << "_phi0." << this->params.phi0;
        return filename.str();
    }

    void setPhi0(double phi0) {
        this->params.phi0 = phi0;
    }
};

#endif //MBL_ED_CAVITYHAMILTONIANGENERATOR_H
