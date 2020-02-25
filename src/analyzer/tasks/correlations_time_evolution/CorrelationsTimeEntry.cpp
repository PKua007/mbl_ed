//
// Created by Piotr Kubala on 21/02/2020.
//

#include "CorrelationsTimeEntry.h"

#include <iterator>

void CorrelationsTimeEntry::Correlations::addObservables(const OccupationEvolution::Occupations &occupations,
                                                         std::size_t borderSize)
{
    std::size_t numberOfSites = occupations.numParticles.size();
    Expects(this->distance < numberOfSites - 2*borderSize);

    double Gsum{};
    for (std::size_t i = borderSize; i < numberOfSites - this->distance - borderSize; i++)
        Gsum += occupations.numParticlesSquared(i, i + this->distance) - occupations.numParticles[i] * occupations.numParticles[i + this->distance];
    std::size_t meanEntries_ = numberOfSites - this->distance - 2 * borderSize;
    this->G += Gsum / meanEntries_;
}

std::string CorrelationsTimeEntry::Correlations::getHeader() const {
    return "G_" + std::to_string(this->distance);
}

std::string CorrelationsTimeEntry::OnsiteFluctuations::getHeader() const {
    return "rho_" + std::to_string(this->i);
}

std::string CorrelationsTimeEntry::Correlations::getValue(std::size_t meanEntries) const {
    std::ostringstream out;
    out << (this->G / meanEntries);
    return out.str();
}

std::string CorrelationsTimeEntry::OnsiteFluctuations::getValue(std::size_t meanEntries) const {
    std::ostringstream out;
    out << (this->rho / meanEntries);
    return out.str();
}

void CorrelationsTimeEntry::OnsiteFluctuations::addObservables(const OccupationEvolution::Occupations &occupations) {
    Expects(this->i < occupations.numParticles.size());

    this->rho += occupations.numParticlesSquared(this->i, this->i) - std::pow(occupations.numParticles[this->i], 2);
}

std::string CorrelationsTimeEntry::CorrelationsTimeEntry::getHeader() const {
    std::ostringstream out;

    out << "t x ";
    auto headerPrinter = [](const auto &entry) { return entry.getHeader(); };
    std::transform(this->correlations.begin(), this->correlations.end(), std::ostream_iterator<std::string>(out, " "),
                   headerPrinter);
    std::transform(this->borderlessCorrelations.begin(), this->borderlessCorrelations.end(),
                   std::ostream_iterator<std::string>(out, " "), headerPrinter);
    std::transform(this->onsiteFluctuations.begin(), this->onsiteFluctuations.end(),
                   std::ostream_iterator<std::string>(out, " "), headerPrinter);

    return out.str();
}

std::string CorrelationsTimeEntry::CorrelationsTimeEntry::getValue() const {
    std::ostringstream out;
    out << this->t << " " << (this->x/this->meanEntries) << " ";

    auto valuePrinter = [this](const auto &entry) { return entry.getValue(meanEntries); };
    std::transform(this->correlations.begin(), this->correlations.end(), std::ostream_iterator<std::string>(out, " "),
                   valuePrinter);
    std::transform(this->borderlessCorrelations.begin(), this->borderlessCorrelations.end(),
                   std::ostream_iterator<std::string>(out, " "), valuePrinter);
    std::transform(this->onsiteFluctuations.begin(), this->onsiteFluctuations.end(),
                   std::ostream_iterator<std::string>(out, " "), valuePrinter);
    return out.str();
}

void CorrelationsTimeEntry::CorrelationsTimeEntry::addObservables(const OccupationEvolution::Occupations &occupations) {
    Expects(this->numberOfSites == occupations.numParticles.size());

    this->meanEntries++;

    for (auto &onsiteFluctuation : this->onsiteFluctuations)
        onsiteFluctuation.addObservables(occupations);
    for (auto &correlation : this->correlations)
        correlation.addObservables(occupations, 0);
    for (auto &borderlessCorrelation : this->borderlessCorrelations)
        borderlessCorrelation.addObservables(occupations, this->borderSize);

    // =, not +=, because this->correlations are already averaged!!! (or rather summed at this point)
    this->x = 2*std::abs(std::accumulate(this->correlations.begin(), this->correlations.end(), 0.,
                                         [](double sum, const Correlations &corr) {
                                             return sum + corr.distance * corr.G;
                                         }));

}

CorrelationsTimeEntry::CorrelationsTimeEntry::CorrelationsTimeEntry(double t, std::size_t borderSize,
                                                                    std::size_t numberOfSites)
        : t{t}, borderSize{borderSize}, numberOfSites{numberOfSites}
{
    Expects(this->numberOfSites > 0);
    Expects(this->numberOfSites - 2*this->borderSize >= 2);

    this->onsiteFluctuations.resize(numberOfSites);
    std::size_t i{};
    for (auto &onsiteFluctuationsEntry : this->onsiteFluctuations)
        onsiteFluctuationsEntry.i = i++;

    this->correlations.resize(numberOfSites - 1);
    std::size_t d = 1;
    for (auto &correlationsEntry : this->correlations)
        correlationsEntry.distance = d++;

    this->borderlessCorrelations.resize(numberOfSites - 1 - 2 * this->borderSize);
    d = 1;
    for (auto &borderlessCorrelationEntry : this->borderlessCorrelations)
        borderlessCorrelationEntry.distance = d++;
}

std::size_t CorrelationsTimeEntry::CorrelationsTimeEntry::getNumberOfSites() const {
    return this->numberOfSites;
}