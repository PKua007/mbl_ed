//
// Created by Piotr Kubala on 21/02/2020.
//

#ifndef MBL_ED_CORRELATIONSTIMEENTRY_H
#define MBL_ED_CORRELATIONSTIMEENTRY_H

#include <vector>

#include "OccupationEvolution.h"

class CorrelationsTimeEntry {
private:
    struct Correlations {
        std::size_t d{};
        double G_d{};

        void addObservables(const OccupationEvolution::Occupations &occupations, std::size_t marginSize_);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string toString(std::size_t numberOfMeanEntries_) const;
    };

    struct OnsiteFluctuations {
        std::size_t i{};
        double rho_i{};

        void addObservables(const OccupationEvolution::Occupations &occupations);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string toString(std::size_t numberOfMeanEntries_) const;
    };

    double t{};
    std::vector<Correlations> correlations{};
    std::vector<Correlations> borderlessCorrelations{};
    std::vector<OnsiteFluctuations> onsiteFluctuations{};

    std::size_t marginSize{};
    std::size_t numberOfSites{};
    std::size_t numberOfMeanEntries{};

    void populateOnsiteFluctuations(std::vector<OnsiteFluctuations> &onsiteFluctuationsVector,
                                    std::size_t numberOfSites_) const;
    void populateCorrelations(std::vector<Correlations> &correlationsVector, size_t numberOfSites_) const;

public:
    CorrelationsTimeEntry() = default;
    CorrelationsTimeEntry(double t, std::size_t marginSize, std::size_t numberOfSites);

    [[nodiscard]] std::size_t getNumberOfSites() const;
    void addObservables(const OccupationEvolution::Occupations &occupations);
    [[nodiscard]] std::string getHeader() const;
    [[nodiscard]] std::string toString() const;
};


#endif //MBL_ED_CORRELATIONSTIMEENTRY_H
