//
// Created by Piotr Kubala on 09/02/2020.
//

#ifndef MBL_ED_ONSITEDISORDER_H
#define MBL_ED_ONSITEDISORDER_H


#include <memory>
#include <algorithm>

#include "utils/Assertions.h"
#include "simulation/DiagonalTerm.h"
#include "simulation/RND.h"

template<typename DisorderGenerator>
class OnsiteDisorder : public DiagonalTerm {
private:
    std::vector<double> onsiteEnergies;
    std::unique_ptr<DisorderGenerator> disorderGenerator;

public:
    explicit OnsiteDisorder(std::unique_ptr<DisorderGenerator> disorderGenerator, std::size_t numberOfSites, RND &rnd)
            : disorderGenerator{std::move(disorderGenerator)}
    {
        Expects(numberOfSites > 0);
        this->onsiteEnergies.resize(numberOfSites);
        this->resampleOnsiteEnergies(rnd);
    }

    /**
     * @brief Clears old onsite disorder values and samples new ones. Used for each new simulation.
     */
    void resampleOnsiteEnergies(RND &rnd) {
        std::generate(this->onsiteEnergies.begin(), this->onsiteEnergies.end(), [this, &rnd]() {
            return (*this->disorderGenerator)(rnd);
        });
    }

    double calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) override {
        Expects(vector.size() == this->onsiteEnergies.size());
        static_cast<void>(generator);

        std::vector<double> elementwiseEnergies;
        elementwiseEnergies.reserve(vector.size());
        std::transform(vector.begin(), vector.end(), this->onsiteEnergies.begin(),
                       std::back_inserter(elementwiseEnergies), std::multiplies<>());
        return std::accumulate(elementwiseEnergies.begin(), elementwiseEnergies.end(), 0., std::plus<>());
    }
};


#endif //MBL_ED_ONSITEDISORDER_H
