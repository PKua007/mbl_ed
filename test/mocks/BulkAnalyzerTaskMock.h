//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_BULKANALYZERTASKMOCK_H
#define MBL_ED_BULKANALYZERTASKMOCK_H

#include <catch2/trompeloeil.hpp>

#include "analyzer/BulkAnalyzerTask.h"


class BulkAnalyzerTaskMock : public trompeloeil::mock_interface<BulkAnalyzerTask> {
public:
    IMPLEMENT_CONST_MOCK1(storeResult);
    IMPLEMENT_MOCK2(analyze);
    IMPLEMENT_CONST_MOCK0(getName);
    IMPLEMENT_CONST_MOCK1(storeState);
    IMPLEMENT_MOCK1(joinRestoredState);
    IMPLEMENT_MOCK0(clear);
};


#endif //MBL_ED_BULKANALYZERTASKMOCK_H
