//
// Created by Piotr Kubala on 21/02/2020.
//

#ifndef MBL_ED_TIMEEVOLUTIONENTRY_H
#define MBL_ED_TIMEEVOLUTIONENTRY_H

#include <utility>
#include <vector>
#include <ostream>

#include "simulation/Restorable.h"
#include "PrimaryObservable.h"
#include "SecondaryObservable.h"


/**
 * @brief Class collecting all stored observable values for a specific moment of time.
 * @details More than one set can be added - it will be avraged.
 */
class TimeEvolutionEntry : public Restorable {
private:
    static constexpr double EPSILON = 1e-10;

    double t{};
    std::vector<double> values{};
    std::size_t numberOfMeanEntries{};

public:
    /**
     * @brief Constructs the class with t = 0 and no values.
     */
    TimeEvolutionEntry() = default;

    /**
     * @brief Constrcts the class with t = @a t and @a numOfValues zero-initialized values (no mean entries yet)
     */
    TimeEvolutionEntry(double t, std::size_t numOfValues) : t{t}, values(numOfValues) { }

    /**
     * @brief Constrcts the class with t = @a t and values initialized from @a values.
     * @details Of course, those @a values are like a single mean entry - they will be properly averaged with other
     * values added with addValues().
     */
    TimeEvolutionEntry(double t, std::vector<double> values)
            : t{t}, values{std::move(values)}, numberOfMeanEntries{1}
    { }

    /**
     * @brief Adds next entry of observable values to be averaged.
     */
    void addValues(const std::vector<double> &newValues);

    friend TimeEvolutionEntry operator+(const TimeEvolutionEntry &first, const TimeEvolutionEntry &second);

    /**
     * @brief Adds all averaging entries from @a other to this instance.
     * @details Times of this instance and other must coincide.
     */
    TimeEvolutionEntry &operator+=(const TimeEvolutionEntry &other);

    friend bool operator==(const TimeEvolutionEntry &first, const TimeEvolutionEntry &second);
    friend std::ostream &operator<<(std::ostream &os, const TimeEvolutionEntry &entry);

    /**
     * @brief Constructs the row of values of all stuff - time and observable values, of course averaged.
     */
    [[nodiscard]] std::string toString() const;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_TIMEEVOLUTIONENTRY_H
