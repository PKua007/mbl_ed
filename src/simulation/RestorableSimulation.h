//
// Created by Piotr Kubala on 31/08/2020.
//

#ifndef MBL_ED_RESTORABLESIMULATION_H
#define MBL_ED_RESTORABLESIMULATION_H

#include <string>

#include "Restorable.h"
#include "SimulationsSpan.h"

class RestorableSimulation : public Restorable {
public:
    virtual void seedRandomGenerators(unsigned long seed) = 0;
    virtual void performSimulation(std::size_t simulationIndex, std::size_t totalSimulations, std::ostream &logger) = 0;
    [[nodiscard]] virtual std::string getTagName() const = 0;
};


#endif //MBL_ED_RESTORABLESIMULATION_H
