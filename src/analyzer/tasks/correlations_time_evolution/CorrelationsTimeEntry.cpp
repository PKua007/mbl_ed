//
// Created by Piotr Kubala on 21/02/2020.
//

#include "CorrelationsTimeEntry.h"

#include <iterator>

void CorrelationsTimeEntry::Correlations::addObservables(const OccupationEvolution::Occupations &occupations,
                                                         std::size_t marginSize_)
{
    std::size_t numberOfSites_ = occupations.numParticles.size();
    Expects(numberOfSites_ > 2 * marginSize_) ;
    std::size_t borderlessNumberOfSites = numberOfSites_ - 2 * marginSize_;
    Expects(this->d < borderlessNumberOfSites);

    double GSiteSum{};
    std::size_t minFirstSite = marginSize_;
    std::size_t maxFirstSite = numberOfSites_ - this->d - marginSize_ - 1;
    for (std::size_t firstSite = minFirstSite; firstSite <= maxFirstSite; firstSite++) {
        std::size_t secondSite = firstSite + this->d;
        double singleG = occupations.numParticlesSquared(firstSite, secondSite) -
                         occupations.numParticles[firstSite] * occupations.numParticles[secondSite];
        GSiteSum += singleG;
    }
    std::size_t meanEntries_ = borderlessNumberOfSites - this->d;
    double GSiteAverage = GSiteSum / meanEntries_;
    this->G_d += GSiteAverage;
}

std::string CorrelationsTimeEntry::Correlations::getHeader() const {
    return "G_" + std::to_string(this->d);
}

std::string CorrelationsTimeEntry::OnsiteFluctuations::getHeader() const {
    return "rho_" + std::to_string(this->i);
}

std::string CorrelationsTimeEntry::Correlations::toString(std::size_t numberOfMeanEntries_) const {
    std::ostringstream out;
    out << (this->G_d / numberOfMeanEntries_);
    return out.str();
}

std::string CorrelationsTimeEntry::OnsiteFluctuations::toString(std::size_t numberOfMeanEntries_) const {
    std::ostringstream out;
    out << (this->rho_i / numberOfMeanEntries_);
    return out.str();
}

void CorrelationsTimeEntry::OnsiteFluctuations::addObservables(const OccupationEvolution::Occupations &occupations) {
    Expects(this->i < occupations.numParticles.size());

    this->rho_i += occupations.numParticlesSquared(this->i, this->i)
                   - std::pow(occupations.numParticles[this->i], 2);
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

std::string CorrelationsTimeEntry::CorrelationsTimeEntry::toString() const {
    std::ostringstream out;
    out << this->t << " " << (this->x/this->numberOfMeanEntries) << " ";

    auto valuePrinter = [this](const auto &entry) { return entry.toString(numberOfMeanEntries); };
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

    this->numberOfMeanEntries++;

    for (auto &onsiteFluctuation : this->onsiteFluctuations)
        onsiteFluctuation.addObservables(occupations);
    for (auto &correlation : this->correlations)
        correlation.addObservables(occupations, 0);
    for (auto &borderlessCorrelation : this->borderlessCorrelations)
        borderlessCorrelation.addObservables(occupations, this->marginSize);

    // =, not +=, because this->correlations are already averaged!!! (or rather summed at this point)
    this->x = 2*std::abs(std::accumulate(this->correlations.begin(), this->correlations.end(), 0.,
                                         [](double sum, const Correlations &corr) {
                                             return sum + corr.d * corr.G_d;
                                         }));

}

CorrelationsTimeEntry::CorrelationsTimeEntry::CorrelationsTimeEntry(double t, std::size_t marginSize,
                                                                    std::size_t numberOfSites)
        : t{t}, marginSize{marginSize}, numberOfSites{numberOfSites}
{
    Expects(this->numberOfSites > 0);
    Expects(this->numberOfSites > 2*this->marginSize);
    std::size_t borderlessNumberOfSites = this->numberOfSites - 2*this->marginSize;
    Expects(borderlessNumberOfSites >= 2);

    this->populateOnsiteFluctuations(this->onsiteFluctuations, this->numberOfSites);
    this->populateCorrelations(this->correlations, this->numberOfSites);
    this->populateCorrelations(this->borderlessCorrelations, borderlessNumberOfSites);
}

void CorrelationsTimeEntry::populateCorrelations(std::vector<Correlations> &correlationsVector,
                                                 std::size_t numberOfSites_) const
{
    correlationsVector.resize(numberOfSites_ - 1);
    std::size_t d = 1;
    for (auto &correlationsEntry : correlationsVector)
        correlationsEntry.d = d++;
}

void CorrelationsTimeEntry::populateOnsiteFluctuations(std::vector<OnsiteFluctuations> &onsiteFluctuationsVector,
                                                       std::size_t numberOfSites_) const
{
    onsiteFluctuationsVector.resize(numberOfSites_);
    std::size_t i{};
    for (auto &onsiteFluctuationsEntry : onsiteFluctuationsVector)
        onsiteFluctuationsEntry.i = i++;
}

std::size_t CorrelationsTimeEntry::CorrelationsTimeEntry::getNumberOfSites() const {
    return this->numberOfSites;
}