//
// Created by Piotr Kubala on 19/02/2020.
//

#ifndef MBL_ED_OCCUPATIONEVOLUTION_H
#define MBL_ED_OCCUPATIONEVOLUTION_H

#include "SymmetricMatrix.h"
#include "simulation/FockBase.h"
#include "simulation/Eigensystem.h"

class OccupationEvolution {
public:
    struct Occupations {
        std::vector<double> numParticles;
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
    [[nodiscard]] static SymmetricMatrix<arma::vec> prepareNumOfParticlesSquaredObservable(const FockBase &fockBase);

    [[nodiscard]] static std::vector<Occupations>
    doPerformEvolution(size_t numSteps, size_t initialFockStateIdx, const FockBase &fockBase,
                       const arma::cx_mat &fockBasisEvolution, const std::vector<arma::vec> &numOfParticlesObservables,
                       const SymmetricMatrix<arma::vec> &numOfParticlesSquaresObservables);

public:
    [[nodiscard]] static std::vector<Occupations> perform(double minTime, double maxTime, std::size_t numSteps,
                                                          std::size_t initialFockStateIdx,
                                                          const Eigensystem &eigensystem);
};


#endif //MBL_ED_OCCUPATIONEVOLUTION_H
