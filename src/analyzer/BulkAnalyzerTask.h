//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_BULKANALYZERTASK_H
#define MBL_ED_BULKANALYZERTASK_H

#include <iosfwd>

#include "AnalyzerTask.h"

/**
 * @brief An AnalyzerTask, which produces long result supposed to be stored in a separate file.
 */
class BulkAnalyzerTask : public AnalyzerTask {
public:
    /**
     * @brief Stores the result of the analyzis task to a given @a out stream.
     */
    virtual void storeResult(std::ostream &out) const = 0;
};

#endif //MBL_ED_BULKANALYZERTASK_H
