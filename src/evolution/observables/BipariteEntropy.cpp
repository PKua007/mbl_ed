//
// Created by Piotr Kubala on 23/09/2020.
//

#include "BipariteEntropy.h"
#include "utils/Assertions.h"
#include "core/FockBasisGenerator.h"

void BipariteEntropy::calculateForState(const arma::cx_vec &state) {
    this->S = 0;

    for (std::size_t particlesOnHalf{}; particlesOnHalf <= this->numOfParticles; particlesOnHalf++) {
        std::size_t particlesOnOtherHalf = this->numOfParticles - particlesOnHalf;
        const auto &halfFockBasis = *this->halfFockBases[particlesOnHalf];
        const auto &otherHalfFockBasis = *this->halfFockBases[particlesOnOtherHalf];
        arma::cx_mat productBasisMatrix(halfFockBasis.size(), otherHalfFockBasis.size());

        for (std::size_t i{}; i < halfFockBasis.size(); i++) {
            for (std::size_t j{}; j < otherHalfFockBasis.size(); j++) {
                FockVector vector = halfFockBasis[i] + otherHalfFockBasis[j];
                std::size_t idx = *(this->fockBasis->findIndex(vector));
                productBasisMatrix(i, j) = state[idx];
            }
        }

        arma::vec singularValues = arma::svd(productBasisMatrix);
        for (double s : singularValues)
            if (s != 0)
                this->S -= s * s * std::log(s * s);
    }
}

BipariteEntropy::BipariteEntropy(std::shared_ptr<FockBasis> fockBasis) : fockBasis{std::move(fockBasis)} {
    std::size_t numOfSites = this->fockBasis->getNumberOfSites();
    Expects(numOfSites % 2 == 0);

    this->numOfParticles = this->fockBasis->getNumberOfParticles();
    std::size_t halfNumOfSites = numOfSites / 2;
    this->halfFockBases.resize(this->numOfParticles + 1);
    for (std::size_t particlesOnHalf{}; particlesOnHalf <= this->numOfParticles; particlesOnHalf++)
        this->halfFockBases[particlesOnHalf] = FockBasisGenerator{}.generate(particlesOnHalf, halfNumOfSites);
}
