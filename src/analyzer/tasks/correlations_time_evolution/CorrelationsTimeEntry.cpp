//
// Created by Piotr Kubala on 21/02/2020.
//

#include "CorrelationsTimeEntry.h"

#include <iterator>

/**
 * @brief The method calculates the site-averaged correlations from given @a occupations and adds the result to G_d.
 * This sum is used in toString() to calculate the average.
 */
void CorrelationsTimeEntry::Correlations::addObservables(const OccupationEvolution::Occupations &occupations)
{
    std::size_t numberOfSites_ = occupations.numParticles.size();
    Expects(numberOfSites_ > 2 * this->marginSize) ;
    std::size_t borderlessNumberOfSites = numberOfSites_ - 2 * this->marginSize;
    Expects(this->d < borderlessNumberOfSites);

    double GSiteSum{};
    std::size_t minFirstSite = this->marginSize;
    std::size_t maxFirstSite = numberOfSites_ - this->d - this->marginSize - 1;
    for (std::size_t firstSite = minFirstSite; firstSite <= maxFirstSite; firstSite++) {
        std::size_t secondSite = firstSite + this->d;
        double singleG = occupations.numParticlesSquared(firstSite, secondSite) -
                         occupations.numParticles[firstSite] * occupations.numParticles[secondSite];
        GSiteSum += singleG;
    }
    std::size_t meanEntries_ = borderlessNumberOfSites - this->d;
    double GSiteAverage = GSiteSum / meanEntries_;
    this->G_d += GSiteAverage;
    this->numberOfMeanEntries++;
}

[[nodiscard]] std::string CorrelationsTimeEntry::Correlations::getHeader() const {
    return "Gm" + std::to_string(this->marginSize) + "_" + std::to_string(this->d);
}

[[nodiscard]] std::string CorrelationsTimeEntry::Correlations::toString() const {
    std::ostringstream out;
    out << (this->G_d / this->numberOfMeanEntries);
    return out.str();
}

/**
 * @brief The method calculates the site-averaged onsite fluctuations from given @a occupations and adds the result to
 * rho_i. This sum is used in toString() to calculate the average.
 */
void CorrelationsTimeEntry::OnsiteFluctuations::addObservables(const OccupationEvolution::Occupations &occupations) {
    Expects(this->i < occupations.numParticles.size());

    this->rho_i += occupations.numParticlesSquared(this->i, this->i)
                   - std::pow(occupations.numParticles[this->i], 2);
    this->numberOfMeanEntries++;
}

std::string CorrelationsTimeEntry::OnsiteFluctuations::getHeader() const {
    return "rho_" + std::to_string(this->i);
}

std::string CorrelationsTimeEntry::OnsiteFluctuations::toString() const {
    std::ostringstream out;
    out << (this->rho_i / this->numberOfMeanEntries);
    return out.str();
}

/**
 * @brief The method calculates the mean onsite occupaitons from given @a occupations and adds the result to n_i. This
 * sum is used in toString() to calculate the average.
 */
void CorrelationsTimeEntry::OnsiteOccupations::addObservables(const OccupationEvolution::Occupations &occupations) {
    Expects(this->i < occupations.numParticles.size());

    this->n_i += occupations.numParticles[this->i];
    this->numberOfMeanEntries++;
}

std::string CorrelationsTimeEntry::OnsiteOccupations::getHeader() const {
    return "n_" + std::to_string(this->i);
}

std::string CorrelationsTimeEntry::OnsiteOccupations::toString() const {
    std::ostringstream out;
    out << (this->n_i / this->numberOfMeanEntries);
    return out.str();
}

void CorrelationsTimeEntry::addObservables(const OccupationEvolution::Occupations &occupations) {
    Expects(this->numberOfSites == occupations.numParticles.size());

    this->numberOfMeanEntries++;

    for (auto &onsiteFluctuation : this->onsiteFluctuations)
        onsiteFluctuation.addObservables(occupations);
    for (auto &onsiteOccupation : this->onsiteOccupations)
        onsiteOccupation.addObservables(occupations);
    for (auto &correlation : this->correlations)
        correlation.addObservables(occupations);
    for (auto &borderlessCorrelation : this->borderlessCorrelations)
        borderlessCorrelation.addObservables(occupations);
}

std::string CorrelationsTimeEntry::getHeader() const {
    std::ostringstream out;

    out << "t ";
    auto headerPrinter = [](const auto &entry) { return entry.getHeader(); };
    std::transform(this->correlations.begin(), this->correlations.end(), std::ostream_iterator<std::string>(out, " "),
                   headerPrinter);
    std::transform(this->borderlessCorrelations.begin(), this->borderlessCorrelations.end(),
                   std::ostream_iterator<std::string>(out, " "), headerPrinter);
    std::transform(this->onsiteFluctuations.begin(), this->onsiteFluctuations.end(),
                   std::ostream_iterator<std::string>(out, " "), headerPrinter);
    std::transform(this->onsiteOccupations.begin(), this->onsiteOccupations.end(),
                   std::ostream_iterator<std::string>(out, " "), headerPrinter);

    return out.str();
}

std::string CorrelationsTimeEntry::toString() const {
    std::ostringstream out;
    out << this->t << " ";

    auto valuePrinter = [](const auto &entry) { return entry.toString(); };
    std::transform(this->correlations.begin(), this->correlations.end(), std::ostream_iterator<std::string>(out, " "),
                   valuePrinter);
    std::transform(this->borderlessCorrelations.begin(), this->borderlessCorrelations.end(),
                   std::ostream_iterator<std::string>(out, " "), valuePrinter);
    std::transform(this->onsiteFluctuations.begin(), this->onsiteFluctuations.end(),
                   std::ostream_iterator<std::string>(out, " "), valuePrinter);
    std::transform(this->onsiteOccupations.begin(), this->onsiteOccupations.end(),
                   std::ostream_iterator<std::string>(out, " "), valuePrinter);
    return out.str();
}

CorrelationsTimeEntry::CorrelationsTimeEntry(double t, std::size_t marginSize, std::size_t numberOfSites)
        : t{t}, numberOfSites{numberOfSites}
{
    Expects(this->numberOfSites > 0);
    Expects(this->numberOfSites > 2*marginSize);
    std::size_t borderlessNumberOfSites = this->numberOfSites - 2*marginSize;
    Expects(borderlessNumberOfSites >= 2);

    this->populateOnsiteFluctuations(this->onsiteFluctuations);
    this->populateOnsiteOccupations(this->onsiteOccupations);
    this->populateCorrelations(this->correlations, 0);
    this->populateCorrelations(this->borderlessCorrelations, marginSize);
}

/**
 * @brief Populates given @a correlationsVector with Correlations for a given margin size.
 */
void CorrelationsTimeEntry::populateCorrelations(std::vector<Correlations> &correlationsVector,
                                                 std::size_t marginSize_) const
{
    std::size_t borderlessNumberOfSites = this->numberOfSites - 2*marginSize_;
    correlationsVector.resize(borderlessNumberOfSites - 1);
    std::size_t d = 1;
    for (auto &correlationsEntry : correlationsVector)
        correlationsEntry = Correlations(d++, marginSize_);
}

/**
 * @brief Populates given @a onsiteFluctuationsVector with OnsiteFluctuations.
 */
void CorrelationsTimeEntry::populateOnsiteFluctuations(std::vector<OnsiteFluctuations> &onsiteFluctuationsVector) const
{
    onsiteFluctuationsVector.resize(this->numberOfSites);
    std::size_t i{};
    for (auto &onsiteFluctuationsEntry : onsiteFluctuationsVector)
        onsiteFluctuationsEntry = OnsiteFluctuations(i++);
}

/**
 * @brief Populates given @a onsiteOccupationsVector with OnsiteOccupations.
 */
void CorrelationsTimeEntry::populateOnsiteOccupations(std::vector<OnsiteOccupations> &onsiteOccupationsVector) const
{
    onsiteOccupationsVector.resize(this->numberOfSites);
    std::size_t i{};
    for (auto &onsiteFluctuationsEntry : onsiteOccupationsVector)
        onsiteFluctuationsEntry = OnsiteOccupations (i++);
}

std::size_t CorrelationsTimeEntry::getNumberOfSites() const {
    return this->numberOfSites;
}