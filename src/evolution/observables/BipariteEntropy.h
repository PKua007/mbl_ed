//
// Created by Piotr Kubala on 23/09/2020.
//

#ifndef MBL_ED_BIPARITEENTROPY_H
#define MBL_ED_BIPARITEENTROPY_H

#include <memory>
#include <vector>

#include "evolution/PrimaryObservable.h"
#include "core/FockBase.h"

class BipariteEntropy : public PrimaryObservable {
private:
    std::size_t numOfSites{};
    std::size_t halfNumOfSites{};
    std::size_t numOfParticles{};
    double S{};
    std::shared_ptr<FockBase> fockBase;
    std::vector<std::unique_ptr<FockBase>> halfFockBases;

public:
    explicit BipariteEntropy(std::shared_ptr<FockBase> fockBase);

    [[nodiscard]] std::vector<std::string> getHeader() const override { return {"S"}; }
    [[nodiscard]] std::vector<double> getValues() const override { return {S}; }
    void calculateForState(const arma::cx_vec &state) override;
};


#endif //MBL_ED_BIPARITEENTROPY_H
