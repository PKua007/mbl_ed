//
// Created by Piotr Kubala on 19/02/2020.
//

#ifndef MBL_ED_OCCUPATIONEVOLUTION_H
#define MBL_ED_OCCUPATIONEVOLUTION_H

#include <valarray>

#include "SymmetricMatrix.h"
#include "simulation/FockBase.h"
#include "Evolver.h"

/**
 * @brief The class performing time evolution of @f$ \left< \hat{n}_i \right> @f$ and
 * @f$ \left< \hat{n}_i \hat{n}_j \right> @f$ expected values of observables, where @f$ \hat{n}_i @f$ is onsite number
 * of parcitles observable.
 */
class OccupationEvolution {
public:
    /**
     * @brief A structure with expected values of @f$ \left< \hat{n}_i \right> @f$ and
     * @f$ \left< \hat{n}_i \hat{n}_j \right> @f$ observables.
     */
    struct Occupations {
        /**
         * @brief Expected values of @f$ \left< \hat{n}_i \right> @f$.
         */
        std::vector<double> numParticles;

        /**
         * @brief Expected values of @f$ \left< \hat{n}_i \hat{n}_j \right> @f$.
         */
        SymmetricMatrix<double> numParticlesSquared;

        Occupations() = default;
        explicit Occupations(std::size_t numberOfSites)
                : numParticles(numberOfSites), numParticlesSquared(numberOfSites)
        { }
    };

private:
    [[nodiscard]] static arma::vec numOfParticlesObservable(const FockBase &fockBase, std::size_t siteIdx);

    [[nodiscard]] static arma::vec numOfParticlesSquaredObservable(const FockBase &fockBase, std::size_t site1Idx,
                                                                   std::size_t site2Idx);

    [[nodiscard]] static double calculateObservableExpectedValue(const arma::vec &observable,
                                                                 const arma::cx_vec &state);

    [[nodiscard]] static std::vector<Occupations> prepareOccupationVector(size_t numSteps, size_t numberOfSites);
    [[nodiscard]] static std::vector<arma::vec> prepareNumOfParticlesObservables(const FockBase &fockBase);

    [[nodiscard]] static SymmetricMatrix<arma::vec>
    prepareNumOfParticlesSquaredObservables(const FockBase &fockBase);

    [[nodiscard]] static std::vector<Occupations>
    doPerformEvolution(std::size_t numSteps, Evolver &evolver, const std::vector<arma::vec> &numOfParticlesObservables,
                       const SymmetricMatrix<arma::vec> &numOfParticlesSquaredObservables, std::ostream &logger);

public:
    /**
     * @brief Perform the evolution of the fock state of index @a initialFockStateIdx, from 0 to @a maxTime, dividing
     * it into @a numSteps steps.
     * @details The evolution operator is created from eigenenergies and eigenstates from @a eigensystem.
     * @return The vector of Occupations, where elements corresponds to expected values of observables in subsequent
     * time steps.
     */
    [[nodiscard]] static std::vector<Occupations> perform(double maxTime, std::size_t numSteps,
                                                          std::size_t initialFockStateIdx,
                                                          const FockBase &fockBase, Evolver &evolver,
                                                          std::ostream &logger);
};


#endif //MBL_ED_OCCUPATIONEVOLUTION_H
