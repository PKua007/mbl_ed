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

/**
 * @brief Parameters used in CavityHamiltonianGenerator. Look there for the description.
 */
struct CavityHamiltonianParameters {
    double J{};
    double U{};
    double U1{};
    double beta{0.5};
    double phi0{};
};

/**
 * @brief Generator of boson hamiltonian from https://arxiv.org/pdf/1902.00357.pdf with modifications.
 *
 * It consist of:
 * <ol>
 *
 * <li> One site tunnelling terms, according to PBC of OBC:
 * \f[ -J \sum_{i=1}^K ( \hat{b}_i^\dagger \hat{b}_{i+1} + \text{c. c.} ) \f]
 *
 * <li> Onsite energy:
 * \f[ \frac{U}{2} \sum_{i=1}^K \hat{n}_i (\hat{n}_i - 1) \f]
 *
 * <li> Onsite disorder potential:
 * \f[ \sum_{i=1}^K E_i \hat{n}_i \f]
 *
 * <li> The long interaction term:
 * \f[ -\frac{U_1}{K} \left( \sum_{i=1}^K \cos(2\pi\beta i+\phi_0) \hat{n}_i \right)^2, \f]
 *
 * </ol>
 * where \f$ i \f$ is the number of site, \f$ K \f$ is the total number of sites, \f$ \hat{b}_i \f$ is annihilation
 * operator, \f$ \hat{n}_i = \hat{b}_i^\dagger\hat{b}_i \f$, \f$ E_i \f$ are random energies sampled by
 * @a DisorderGenerator and the rest of the letters are parameters from CavityHamiltonianParameters. If
 * \f$ \beta=0.5 \f$, \f$ \phi_0=0 \f$, then the hamiltonian has the same form as in the reference above.
 *
 * @tparam DisorderGenerator class whose operator() returns random energies as double
 */
template<typename DisorderGenerator>
class CavityHamiltonianGenerator : public HamiltonianGenerator {
private:
    std::vector<double> onsiteEnergies;
    std::unique_ptr<DisorderGenerator> disorderGenerator;
    CavityHamiltonianParameters params;

    /* All diagonal terms below are described above the class */

    [[nodiscard]] double getOnsiteDisorderEnergy(const FockBase::Vector &vector) const {
        std::vector<double> elementwiseEnergies;
        elementwiseEnergies.reserve(vector.size());
        std::transform(vector.begin(), vector.end(), this->onsiteEnergies.begin(),
                       std::back_inserter(elementwiseEnergies), std::multiplies<>());
        return std::accumulate(elementwiseEnergies.begin(), elementwiseEnergies.end(), 0., std::plus<>());
    }

    [[nodiscard]] double getOnsiteEnergy(const FockBase::Vector &vector) const {
        auto bosonAccumulator = [](auto sum, auto numberOfParticles) {
            return sum + numberOfParticles*(numberOfParticles - 1);
        };
        return this->params.U / 2 * std::accumulate(vector.begin(), vector.end(), 0., bosonAccumulator);
    }

    [[nodiscard]] double getLongInteractionEnergy(const FockBase::Vector &vector) const {
        std::size_t elementIndex{};
        auto plusMinusAccumulator = [&elementIndex, this](auto sum, auto element) {
            return sum + std::cos(2*M_PI*params.beta*(elementIndex++) + params.phi0) * element;
        };
        double populationImbalance = std::accumulate(vector.begin(), vector.end(), 0., plusMinusAccumulator);
        return -this->params.U1 / vector.size() * populationImbalance * populationImbalance;
    }

public:
    /**
     * @brief Creates the hamiltonian generator. Onsite disorder energies are initialized with random generator.
     */
    CavityHamiltonianGenerator(std::unique_ptr<FockBase> fockBase, CavityHamiltonianParameters params,
                               std::unique_ptr<DisorderGenerator> disorderGenerator, bool usePbc = true)
            : HamiltonianGenerator(std::move(fockBase), usePbc), disorderGenerator(std::move(disorderGenerator)),
              params{params}
    {
        this->resampleOnsiteEnergies();
    }

    /**
     * @brief Clears old onsite disorder values and samples new ones. Used for each new simulation.
     */
    void resampleOnsiteEnergies() {
        this->onsiteEnergies.resize(this->fockBase->getNumberOfSites());
        std::generate(this->onsiteEnergies.begin(), this->onsiteEnergies.end(), std::ref(*(this->disorderGenerator)));
    }

    [[nodiscard]] const std::vector<double> &getOnsiteEnergies() const {
        return this->onsiteEnergies;
    }

    [[nodiscard]] double getDiagonalElement(const FockBase::Vector &vector) const override {
        return getOnsiteDisorderEnergy(vector) + getOnsiteEnergy(vector) + getLongInteractionEnergy(vector);
    }

    /**
     * @brief Returns constant -J (see CavityHamiltonianParameters) for a hop between neighbouring
     * sites (according to PBC or OBC) and 0 for larger hops.
     */
    [[nodiscard]] double getHoppingTerm(std::size_t fromSiteIndex, std::size_t toSiteIndex) const override {
        Expects(this->getSiteDistance(fromSiteIndex, toSiteIndex) == 1);

        return -this->params.J;
    }

    void setPhi0(double phi0) {
        this->params.phi0 = phi0;
    }
};

#endif //MBL_ED_CAVITYHAMILTONIANGENERATOR_H
