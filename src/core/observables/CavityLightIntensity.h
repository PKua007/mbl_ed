//
// Created by Piotr Kubala on 06/12/2020.
//

#ifndef MBL_ED_CAVITYLIGHTINTENSITY_H
#define MBL_ED_CAVITYLIGHTINTENSITY_H


#include "core/SecondaryObservable.h"

class CavityLightIntensity : public SecondaryObservable {
private:
    double lightIntensity{};

public:
    [[nodiscard]] std::vector<std::string> getHeader() const override;
    [[nodiscard]] std::vector<double> getValues() const override;
    void calculateForObservables(const std::vector<std::shared_ptr<PrimaryObservable>> &primaryObservables) override;
};


#endif //MBL_ED_CAVITYLIGHTINTENSITY_H
