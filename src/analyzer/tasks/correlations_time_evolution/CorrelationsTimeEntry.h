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
        std::size_t distance{};
        double G{};

        void addObservables(const OccupationEvolution::Observables &observables, std::size_t borderSize);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string getValue(std::size_t meanEntries) const;
    };

    struct OnsiteFluctuations {
        std::size_t i{};
        double rho{};

        void addObservables(const OccupationEvolution::Observables &observables);
        [[nodiscard]] std::string getHeader() const;
        [[nodiscard]] std::string getValue(std::size_t meanEntries) const;
    };

    double t{};
    double x{};
    std::vector<Correlations> correlations{};
    std::vector<Correlations> borderlessCorrelations{};
    std::vector<OnsiteFluctuations> onsiteFluctuations{};

    std::size_t borderSize{};
    std::size_t numberOfSites{};
    std::size_t meanEntries{};

public:
    CorrelationsTimeEntry() = default;
    CorrelationsTimeEntry(double t, std::size_t borderSize, std::size_t numberOfSites);

    [[nodiscard]] std::size_t getNumberOfSites() const;
    void addObservables(const OccupationEvolution::Observables &observables);
    [[nodiscard]] std::string getHeader() const;
    [[nodiscard]] std::string getValue() const;
};


#endif //MBL_ED_CORRELATIONSTIMEENTRY_H
