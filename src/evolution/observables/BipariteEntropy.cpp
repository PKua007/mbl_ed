//
// Created by Piotr Kubala on 23/09/2020.
//

#include "BipariteEntropy.h"
#include "utils/Assertions.h"
#include "core/FockBaseGenerator.h"

void BipariteEntropy::calculateForState(const arma::cx_vec &state) {
    this->S = 0;

    for (std::size_t particlesOnHalf{}; particlesOnHalf <= this->numOfParticles; particlesOnHalf++) {
        std::size_t particlesOnOtherHalf = this->numOfParticles - particlesOnHalf;
        const auto &halfFockBase = *this->halfFockBases[particlesOnHalf];
        const auto &otherHalfFockBase = *this->halfFockBases[particlesOnOtherHalf];
        arma::cx_mat productBasisMatrix(halfFockBase.size(), otherHalfFockBase.size());

        for (std::size_t i{}; i < halfFockBase.size(); i++) {
            for (std::size_t j{}; j < otherHalfFockBase.size(); j++) {
                FockVector vector = halfFockBase[i] + otherHalfFockBase[j];
                std::size_t idx = *(this->fockBase->findIndex(vector));
                productBasisMatrix(i, j) = state[idx];
            }
        }

        arma::vec singularValues = arma::svd(productBasisMatrix);
        for (double s : singularValues)
            if (s != 0)
                this->S -= s * s * std::log(s * s);
    }
}

BipariteEntropy::BipariteEntropy(std::shared_ptr<FockBase> fockBase) : fockBase{std::move(fockBase)} {
    std::size_t numOfSites = this->fockBase->getNumberOfSites();
    Expects(numOfSites % 2 == 0);

    this->numOfParticles = this->fockBase->getNumberOfParticles();
    std::size_t halfNumOfSites = numOfSites / 2;
    this->halfFockBases.resize(this->numOfParticles + 1);
    for (std::size_t particlesOnHalf{}; particlesOnHalf <= this->numOfParticles; particlesOnHalf++)
        this->halfFockBases[particlesOnHalf] = FockBaseGenerator{}.generate(particlesOnHalf, halfNumOfSites);
}