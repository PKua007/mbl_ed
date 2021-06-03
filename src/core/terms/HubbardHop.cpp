//
// Created by Piotr Kubala on 09/02/2020.
//

#include <ZipIterator.hpp>

#include "HubbardHop.h"
#include "core/HamiltonianGenerator.h"

double HubbardHop::calculate(const HopData &hopData, const HamiltonianGenerator &generator) const {
    std::size_t distance = generator.getSiteDistance(hopData.fromSite, hopData.toSite);
    std::size_t paramIdx = std::find(this->hoppingDistances.begin(), this->hoppingDistances.end(), distance)
                           - this->hoppingDistances.begin();
    Expects(paramIdx != this->hoppingDistances.size());

    return -this->Js[paramIdx];
}

HubbardHop::HubbardHop(double J) : Js{J}, hoppingDistances{1} {
    Expects(J != 0);
}

HubbardHop::HubbardHop(std::vector<double> Js) : Js{std::move(Js)} {
    Expects(!this->Js.empty());
    this->hoppingDistances.resize(this->Js.size());
    std::iota(this->hoppingDistances.begin(), this->hoppingDistances.end(), 1);
}

HubbardHop::HubbardHop(std::vector<std::size_t> hoppingDistances, std::vector<double> Js)
        : Js{std::move(Js)}, hoppingDistances{std::move(hoppingDistances)}
{
    Expects(!this->Js.empty());
    Expects(this->Js.size() == this->hoppingDistances.size());
    Expects(std::all_of(this->hoppingDistances.begin(), this->hoppingDistances.end(),
                        [](std::size_t dist) { return dist > 0; }));

    auto zipped = Zip(this->hoppingDistances, this->Js);
    std::sort(zipped.begin(), zipped.end());
    Expects(std::unique(this->hoppingDistances.begin(), this->hoppingDistances.end()) == this->hoppingDistances.end());
}
