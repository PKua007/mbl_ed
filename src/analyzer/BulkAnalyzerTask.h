//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_BULKANALYZERTASK_H
#define MBL_ED_BULKANALYZERTASK_H

#include <iosfwd>

#include "AnalyzerTask.h"

class BulkAnalyzerTask : public AnalyzerTask {
public:
    virtual void storeResult(std::ostream &out) const = 0;
};

#endif //MBL_ED_BULKANALYZERTASK_H
