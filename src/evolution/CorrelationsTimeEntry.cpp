//
// Created by Piotr Kubala on 21/02/2020.
//

#include "CorrelationsTimeEntry.h"

#include "simulation/RestorableHelper.h"

#include <iterator>


void CorrelationsTimeEntry::addValues(const std::vector<double> &newValues) {
    Expects(this->values.size() == newValues.size());

    std::transform(newValues.begin(), newValues.end(), this->values.begin(), this->values.begin(), std::plus{});
    this->numberOfMeanEntries++;
}

std::string CorrelationsTimeEntry::toString() const {
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

void CorrelationsTimeEntry::storeState(std::ostream &binaryOut) const {
    binaryOut.write(reinterpret_cast<const char*>(&this->t), sizeof(this->t));
    binaryOut.write(reinterpret_cast<const char*>(&this->numberOfMeanEntries), sizeof(this->numberOfMeanEntries));
    Assert(binaryOut.good());
    RestorableHelper::storeStateForVector(this->values, binaryOut);
}

void CorrelationsTimeEntry::joinRestoredState(std::istream &binaryIn) {
    double tRestored{};
    std::size_t numberOfMeanEntriesRestored{};
    binaryIn.read(reinterpret_cast<char*>(&tRestored), sizeof(tRestored));
    binaryIn.read(reinterpret_cast<char*>(&numberOfMeanEntriesRestored), sizeof(numberOfMeanEntriesRestored));
    Assert(binaryIn.good());

    CorrelationsTimeEntry restoredEntry;
    // restoredEntry has an empty vector at the beginning, so it will correctly populate it
    RestorableHelper::joinRestoredStateForVector(restoredEntry.values, binaryIn);
    restoredEntry.t = tRestored;
    restoredEntry.numberOfMeanEntries = numberOfMeanEntriesRestored;

    (*this) += restoredEntry;
}

void CorrelationsTimeEntry::clear() {
    this->numberOfMeanEntries = 0;
    std::fill(this->values.begin(), this->values.end(), 0.);
}

CorrelationsTimeEntry &CorrelationsTimeEntry::operator+=(const CorrelationsTimeEntry &other) {
    Expects(this->values.size() == other.values.size());
    Expects(std::abs(this->t - other.t) < EPSILON);

    std::transform(other.values.begin(), other.values.end(), this->values.begin(), this->values.begin(), std::plus{});
    this->numberOfMeanEntries += other.numberOfMeanEntries;
    return *this;
}

CorrelationsTimeEntry operator+(const CorrelationsTimeEntry &first, const CorrelationsTimeEntry &second) {
    CorrelationsTimeEntry result(first);
    result += second;
    return result;
}

bool operator==(const CorrelationsTimeEntry &first, const CorrelationsTimeEntry &second) {
    return first.values == second.values && std::abs(first.t - second.t) < CorrelationsTimeEntry::EPSILON;
}

std::ostream &operator<<(std::ostream &out, const CorrelationsTimeEntry &entry) {
    return out << entry.toString();
}
