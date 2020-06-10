//
// Created by pkua on 08.06.2020.
//

#include "OnsiteDisorder.h"

OnsiteDisorder::OnsiteDisorder(std::unique_ptr<DisorderGenerator> disorderGenerator, std::size_t numberOfSites,
                               RND &rnd)
        : disorderGenerator{std::move(disorderGenerator)}
{
    Expects(numberOfSites > 0);
    this->onsiteEnergies.resize(numberOfSites);
    this->resampleOnsiteEnergies(rnd);
}

void OnsiteDisorder::resampleOnsiteEnergies(RND &rnd) {
    std::generate(this->onsiteEnergies.begin(), this->onsiteEnergies.end(), [this, &rnd]() {
        return this->disorderGenerator->generate(rnd);
    });
}

double OnsiteDisorder::calculate(const FockBase::Vector &vector, const HamiltonianGenerator &generator) {
    Expects(vector.size() == this->onsiteEnergies.size());
    static_cast<void>(generator);

    std::vector<double> elementwiseEnergies;
    elementwiseEnergies.reserve(vector.size());
    std::transform(vector.begin(), vector.end(), this->onsiteEnergies.begin(),
                   std::back_inserter(elementwiseEnergies), std::multiplies<>());
    return std::accumulate(elementwiseEnergies.begin(), elementwiseEnergies.end(), 0., std::plus<>());
}
