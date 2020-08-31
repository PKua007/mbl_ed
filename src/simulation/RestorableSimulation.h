//
// Created by Piotr Kubala on 31/08/2020.
//

#ifndef MBL_ED_RESTORABLESIMULATION_H
#define MBL_ED_RESTORABLESIMULATION_H

#include "Restorable.h"
#include "SimulationsSpan.h"

class RestorableSimulation : public Restorable {
public:
    virtual void seedRandomGenerators(unsigned long seed) = 0;
    virtual void performSimulations(const SimulationsSpan &simulationsSpan, std::ostream &logger) = 0;
};


#endif //MBL_ED_RESTORABLESIMULATION_H
