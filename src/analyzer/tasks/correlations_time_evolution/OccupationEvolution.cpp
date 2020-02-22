//
// Created by Piotr Kubala on 19/02/2020.
//

#include "OccupationEvolution.h"

std::vector<OccupationEvolution::Observables> OccupationEvolution::perform(const std::vector<double> &times,
                                                                           size_t initialIdx,
                                                                           const Eigensystem &eigensystem)
{
    Expects(eigensystem.hasFockBase());
#ifndef NDEBUG
    Expects(eigensystem.isOrthonormal());
#endif

    const auto &fockBase = eigensystem.getFockBase();
    const auto &eigenvectors = eigensystem.getEigenstates();
    const auto &eigenenergies = eigensystem.getEigenenergies();
    std::size_t numberOfSites = fockBase.getNumberOfSites();

    std::vector<Observables> observablesEvolution;
    observablesEvolution.resize(times.size());
    std::transform(times.begin(), times.end(), observablesEvolution.begin(),
                   [numberOfSites](double time) { return Observables(numberOfSites, time); });

    for (std::size_t site{}; site < numberOfSites; site++) {
        arma::mat matrixElements = numberOfParticlesObservable(fockBase, eigenvectors, site);
        arma::mat evolutionTerms = calculateEvolutionTerms(std::move(matrixElements), fockBase,
                                                                 eigenvectors, initialIdx);

        for (std::size_t timeIdx{}; timeIdx < times.size(); timeIdx++) {
            double time = times[timeIdx];
            double value = calculateObservableValue(evolutionTerms, time, eigenenergies);
            observablesEvolution[timeIdx].ns[site] = value;
        }
    }

    for (std::size_t site1{}; site1 < numberOfSites; site1++) {
        for (std::size_t site2 = site1; site2 < numberOfSites; site2++) {
            arma::mat matrixElements = numberOfParticlesSquaredObservable(fockBase, eigenvectors,
                                                                                site1, site2);
            arma::mat evolutionTerms = calculateEvolutionTerms(std::move(matrixElements), fockBase,
                                                                     eigenvectors, initialIdx);

            for (std::size_t timeIdx{}; timeIdx < times.size(); timeIdx++) {
                double time = times[timeIdx];
                double value = calculateObservableValue(evolutionTerms, time, eigenenergies);
                observablesEvolution[timeIdx].nns(site1, site2) = value;
            }
        }
    }

    return observablesEvolution;
}

double OccupationEvolution::calculateObservableValue(const arma::mat &evolutionTerms, double time,
                                                     const arma::vec &eigenenergies)
{
    double value{};

    for (std::size_t elemI{}; elemI < evolutionTerms.n_rows; elemI++) {
        value += evolutionTerms(elemI, elemI);
        for (std::size_t elemJ = elemI + 1; elemJ < evolutionTerms.n_rows; elemJ++)
            value += 2 * std::cos((eigenenergies[elemI] - eigenenergies[elemJ]) * time) * evolutionTerms(elemI, elemJ);
    }
    return value;
}

arma::mat OccupationEvolution::calculateEvolutionTerms(arma::mat matrixElements, const FockBase &fockBase,
                                                             const arma::mat &eigenvectors, std::size_t initialIdx)
{
    for (std::size_t elemI{}; elemI < fockBase.size(); elemI++)
        for (std::size_t elemJ = elemI; elemJ < fockBase.size(); elemJ++)
            matrixElements(elemI, elemJ) *= eigenvectors(initialIdx, elemI) * eigenvectors(initialIdx, elemJ);
    return matrixElements;
}

arma::mat OccupationEvolution::numberOfParticlesObservable(const FockBase &fockBase,
                                                                 const arma::mat &eigenvectors, std::size_t site)
{
//    arma::mat matrixElements(fockBase.size(), fockBase.size(), arma::fill::zeros);
//    for (std::size_t elemI{}; elemI < fockBase.size(); elemI++)
//        for (std::size_t elemJ = elemI; elemJ < fockBase.size(); elemJ++)
//            for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++)
//                matrixElements(elemI, elemJ) += eigenvectors(fockIdx, elemI) * eigenvectors(fockIdx, elemJ) * fockBase[fockIdx][site];
    arma::vec vec(fockBase.size());
    for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++)
        vec[fockIdx] = fockBase[fockIdx][site];

    return eigenvectors.t() * arma::diagmat(vec) * eigenvectors;
}

arma::mat OccupationEvolution::numberOfParticlesSquaredObservable(const FockBase &fockBase,
                                                                        const arma::mat &eigenvectors,
                                                                        std::size_t site1, std::size_t site2)
{
//    arma::mat matrixElements(fockBase.size(), fockBase.size(), arma::fill::zeros);
//    for (std::size_t elemI{}; elemI < fockBase.size(); elemI++) {
//        for (std::size_t elemJ = elemI; elemJ < fockBase.size(); elemJ++) {
//            for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++) {
//                matrixElements(elemI, elemJ) += eigenvectors(fockIdx, elemI) * eigenvectors(fockIdx, elemJ)
//                                                * fockBase[fockIdx][site1] * fockBase[fockIdx][site2];
//            }
//        }
//    }
//    return matrixElements;
    arma::vec vec(fockBase.size());
    for (std::size_t fockIdx{}; fockIdx < fockBase.size(); fockIdx++)
        vec[fockIdx] = fockBase[fockIdx][site1] * fockBase[fockIdx][site2];

    return eigenvectors.t() * arma::diagmat(vec) * eigenvectors;
}