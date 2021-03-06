//
// Created by Piotr Kubala on 23/09/2020.
//

#ifndef MBL_ED_BIPARITEENTROPY_H
#define MBL_ED_BIPARITEENTROPY_H

#include <memory>
#include <vector>

#include "core/PrimaryObservable.h"
#include "core/FockBasis.h"

/**
 * @brief The observable defined as -Tr [rho_a log(rho_a)].
 * @details rho_A is the density matrix of a whole system traced over the "right" half of sites, so it is the density
 * matrix of the "left" half of sites. The class has only one field - this entropy. If the number of sites is odd,
 * "right" part is (n-1)/2, and "left" is (n+1)/2.
 */
class BipariteEntropy : public PrimaryObservable {
private:
    std::size_t numOfParticles{};
    double S{};
    std::shared_ptr<FockBasis> fockBasis;
    std::vector<std::shared_ptr<FockBasis>> halfFockBases;
    std::vector<std::shared_ptr<FockBasis>> otherHalfFockBases;

public:
    explicit BipariteEntropy(std::shared_ptr<FockBasis> fockBasis);

    [[nodiscard]] std::vector<std::string> getHeader() const override { return {"S"}; }
    [[nodiscard]] std::vector<double> getValues() const override { return {S}; }
    void calculateForState(const arma::cx_vec &state) override;
};


#endif //MBL_ED_BIPARITEENTROPY_H
