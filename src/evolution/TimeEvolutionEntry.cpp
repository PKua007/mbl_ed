//
// Created by Piotr Kubala on 21/02/2020.
//

#include "TimeEvolutionEntry.h"

#include "simulation/RestorableHelper.h"

#include <iterator>


void TimeEvolutionEntry::addValues(const std::vector<double> &newValues) {
    Expects(this->values.size() == newValues.size());

    std::transform(newValues.begin(), newValues.end(), this->values.begin(), this->values.begin(), std::plus{});
    this->numberOfMeanEntries++;
}

std::string TimeEvolutionEntry::toString() const {
    std::ostringstream out;
    out << this->t << " ";
    if (this->numberOfMeanEntries != 0) {
        std::transform(this->values.begin(), this->values.end(), std::ostream_iterator<double>(out, " "),
                       [this](double value) { return value / this->numberOfMeanEntries; });
    } else {
        std::copy(this->values.begin(), this->values.end(), std::ostream_iterator<double>(out, " "));
    }
    return out.str();
}

void TimeEvolutionEntry::storeState(std::ostream &binaryOut) const {
    binaryOut.write(reinterpret_cast<const char*>(&this->t), sizeof(this->t));
    binaryOut.write(reinterpret_cast<const char*>(&this->numberOfMeanEntries), sizeof(this->numberOfMeanEntries));
    Assert(binaryOut.good());
    RestorableHelper::storeStateForVector(this->values, binaryOut);
}

void TimeEvolutionEntry::joinRestoredState(std::istream &binaryIn) {
    double tRestored{};
    std::size_t numberOfMeanEntriesRestored{};
    binaryIn.read(reinterpret_cast<char*>(&tRestored), sizeof(tRestored));
    binaryIn.read(reinterpret_cast<char*>(&numberOfMeanEntriesRestored), sizeof(numberOfMeanEntriesRestored));
    Assert(binaryIn.good());

    TimeEvolutionEntry restoredEntry;
    // restoredEntry has an empty vector at the beginning, so it will correctly populate it
    RestorableHelper::joinRestoredStateForVector(restoredEntry.values, binaryIn);
    restoredEntry.t = tRestored;
    restoredEntry.numberOfMeanEntries = numberOfMeanEntriesRestored;

    (*this) += restoredEntry;
}

void TimeEvolutionEntry::clear() {
    this->numberOfMeanEntries = 0;
    std::fill(this->values.begin(), this->values.end(), 0.);
}

TimeEvolutionEntry &TimeEvolutionEntry::operator+=(const TimeEvolutionEntry &other) {
    Expects(this->values.size() == other.values.size());
    Expects(std::abs(this->t - other.t) < EPSILON);

    std::transform(other.values.begin(), other.values.end(), this->values.begin(), this->values.begin(), std::plus{});
    this->numberOfMeanEntries += other.numberOfMeanEntries;
    return *this;
}

TimeEvolutionEntry operator+(const TimeEvolutionEntry &first, const TimeEvolutionEntry &second) {
    TimeEvolutionEntry result(first);
    result += second;
    return result;
}

bool operator==(const TimeEvolutionEntry &first, const TimeEvolutionEntry &second) {
    return first.values == second.values && std::abs(first.t - second.t) < TimeEvolutionEntry::EPSILON;
}

std::ostream &operator<<(std::ostream &out, const TimeEvolutionEntry &entry) {
    return out << entry.toString();
}
