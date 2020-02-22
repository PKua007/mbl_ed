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
    struct Observables {
        double time{};
        std::vector<double> ns;
        SymmetricMatrix nns;

        Observables() = default;
        explicit Observables(std::size_t numberOfSites, double time)
                : time{time}, ns(numberOfSites), nns(numberOfSites)
        { }
    };

private:
    [[nodiscard]] static arma::mat numberOfParticlesObservable(const FockBase &fockBase,
                                                                     const arma::mat &eigenvectors, std::size_t site);

    [[nodiscard]] static arma::mat numberOfParticlesSquaredObservable(const FockBase &fockBase,
                                                                            const arma::mat &eigenvectors,
                                                                            std::size_t site1, std::size_t site2);

    [[nodiscard]] static arma::mat calculateEvolutionTerms(arma::mat matrixElements,
                                                                 const FockBase &fockBase,
                                                                 const arma::mat &eigenvectors, std::size_t initialIdx);

    [[nodiscard]] static double calculateObservableValue(const arma::mat &evolutionTerms, double time,
                                                         const arma::vec &eigenenergies);

public:
    [[nodiscard]] static std::vector<Observables> perform(const std::vector<double> &times, size_t initialIdx,
                                                          const Eigensystem &eigensystem);
};


#endif //MBL_ED_OCCUPATIONEVOLUTION_H
