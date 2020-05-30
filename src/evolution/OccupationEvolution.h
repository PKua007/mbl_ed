//
// Created by Piotr Kubala on 19/02/2020.
//

#ifndef MBL_ED_OCCUPATIONEVOLUTION_H
#define MBL_ED_OCCUPATIONEVOLUTION_H

#include <memory>
#include <utility>

#include "SymmetricMatrix.h"
#include "simulation/FockBase.h"
#include "Evolver.h"
#include "EvolutionTimeSegment.h"

/**
 * @brief The class performing time evolution of @f$ \left< \hat{n}_i \right> @f$ and
 * @f$ \left< \hat{n}_i \hat{n}_j \right> @f$ expected values of observables, where @f$ \hat{n}_i @f$ is onsite number
 * of particles observable.
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
    std::shared_ptr<FockBase> fockBase;
    std::vector<arma::vec> numOfParticlesObservables;
    SymmetricMatrix<arma::vec> numOfParticlesSquaredObservables;
    std::size_t timeStep{};

    [[nodiscard]] arma::vec calculateNumOfParticlesObservable(std::size_t siteIdx) const;
    [[nodiscard]] arma::vec calculateNumOfParticlesSquaredObservable(std::size_t site1Idx, std::size_t site2Idx) const;
    [[nodiscard]] double calculateObservableExpectedValue(const arma::vec &observable, const arma::cx_vec &state) const;
    [[nodiscard]] std::vector<Occupations> prepareOccupationVector(size_t numSteps, size_t numberOfSites);
    void prepareNumOfParticlesObservables();
    void prepareNumOfParticlesSquaredObservables();
    [[nodiscard]] std::vector<Occupations> performTimeSegmentEvolution(std::size_t numSteps, Evolver &evolver,
                                                                       std::ostream &logger);

public:
    explicit OccupationEvolution(std::shared_ptr<FockBase> fockBase);

    /**
     * @brief Perform the evolution of the fock state of index @a initialFockStateIdx for times specified by
     * @a timeSegmentation.
     * @details The actual evolution is performed by a given Evolver. Time segmentation is described in
     * CorrelationsTimeEvolutionParameters.
     * @return The vector of Occupations, where elements corresponds to expected values of observables in subsequent
     * time steps.
     */
    [[nodiscard]] std::vector<Occupations> perform(const std::vector<EvolutionTimeSegment> &timeSegmentation,
                                                   std::size_t initialFockStateIdx, Evolver &evolver,
                                                   std::ostream &logger);
};


#endif //MBL_ED_OCCUPATIONEVOLUTION_H
