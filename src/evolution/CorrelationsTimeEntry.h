//
// Created by Piotr Kubala on 21/02/2020.
//

#ifndef MBL_ED_CORRELATIONSTIMEENTRY_H
#define MBL_ED_CORRELATIONSTIMEENTRY_H

#include <utility>
#include <vector>
#include <ostream>

#include "simulation/Restorable.h"
#include "PrimaryObservable.h"
#include "SecondaryObservable.h"


/**
 * @brief Class collecting OccupationEvolution::Occupations for subsequent eigensystems and calculating disorder
 * averages of correlations as well as onsite occupations and fluctuations calculated from before mentioned occupations
 * for a specific moment of time.
 * @details The class contains:
 * <ul>
 *     <li>site-averaged correlations, for all distanced from 1 to <em> number of sites - 1 </em>
 *     (see CorrelationsTimeEntry::Correlations in the code)
 *     <li>site-averaget correlations, but taking into account only sites after removing a few on the borders
 *     (@a marginSize parameters from the constructor)
 *     <li>onsite fluctuations (see CorrelationsTimeEntry::OnsiteFluctuations in the code)
 *     <li>onsite occupaitons (see CorrelationsTimeEntry::OnsiteOccupations in the code)
 * </ul>
 */
class CorrelationsTimeEntry : public Restorable {
private:
    static constexpr double EPSILON = 1e-10;

    double t{};
    std::vector<double> values{};
    std::size_t numberOfMeanEntries{};

public:
    CorrelationsTimeEntry() = default;
    CorrelationsTimeEntry(double t, std::size_t numOfValues) : t{t}, values(numOfValues) { }
    CorrelationsTimeEntry(double t, std::vector<double> values)
            : t{t}, values{std::move(values)}, numberOfMeanEntries{1}
    { }

    /**
     * @brief Adds next entry of observables
     */
    void addValues(const std::vector<double> &newValues);

    friend CorrelationsTimeEntry operator+(const CorrelationsTimeEntry &first, const CorrelationsTimeEntry &second);
    CorrelationsTimeEntry &operator+=(const CorrelationsTimeEntry &other);
    friend bool operator==(const CorrelationsTimeEntry &first, const CorrelationsTimeEntry &second);
    friend std::ostream &operator<<(std::ostream &os, const CorrelationsTimeEntry &entry);

    /**
     * @brief Constructs the row of values of all stuff - time and observable values
     */
    [[nodiscard]] std::string toString() const;

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;
};


#endif //MBL_ED_CORRELATIONSTIMEENTRY_H
