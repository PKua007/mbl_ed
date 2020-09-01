//
// Created by Piotr Kubala on 31/08/2020.
//

#ifndef MBL_ED_RESTORABLESIMULATIONEXECUTOR_H
#define MBL_ED_RESTORABLESIMULATIONEXECUTOR_H

#include <memory>

#include "RestorableSimulation.h"

class RestorableSimulationExecutor {
private:
    SimulationsSpan simulationsSpan;
    std::string fileSignature;
    bool workloadSplit{};
    bool shouldSave_{};

    void storeSimulations(const RestorableSimulation &simulation, std::ostream &binaryOut,
                          std::size_t simulationIndex, bool finished) const;
    std::pair<bool, std::size_t> restoreSimulations(RestorableSimulation &simulation, std::istream &binaryIn) const;

public:
    explicit RestorableSimulationExecutor(const SimulationsSpan &simulationsSpan,
                                          std::string fileSignature, bool workloadSplit);

    void performSimulations(RestorableSimulation &simulation, unsigned long seed, std::ostream &logger);
    [[nodiscard]] bool shouldSave() const { return this->shouldSave_; }
};


#endif //MBL_ED_RESTORABLESIMULATIONEXECUTOR_H
