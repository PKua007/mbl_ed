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
        std::vector<double> ns;
        SymmetricMatrix nns;

        Observables() = default;
        explicit Observables(std::size_t numberOfSites)
                : ns(numberOfSites), nns(numberOfSites)
        { }
    };

private:
    [[nodiscard]] static arma::vec numberOfParticlesObservable(const FockBase &fockBase,
                                                                     const arma::mat &eigenvectors, std::size_t site);

    [[nodiscard]] static arma::vec numberOfParticlesSquaredObservable(const FockBase &fockBase,
                                                                            const arma::mat &eigenvectors,
                                                                            std::size_t site1, std::size_t site2);

    [[nodiscard]] static double calculateObservableValue(const arma::vec &observable, const arma::cx_vec &state);

public:
    [[nodiscard]] static std::vector<Observables> perform(double minTime, double maxTime, std::size_t numSteps,
                                                        size_t initialIdx,
                                                          const Eigensystem &eigensystem);
};


#endif //MBL_ED_OCCUPATIONEVOLUTION_H
