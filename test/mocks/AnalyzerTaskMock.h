//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_ANALYZERTASKMOCK_H
#define MBL_ED_ANALYZERTASKMOCK_H

#include <catch2/trompeloeil.hpp>

#include "analyzer/AnalyzerTask.h"

class AnalyzerTaskMock : public trompeloeil::mock_interface<AnalyzerTask> {
public:
    IMPLEMENT_MOCK2(analyze);
    IMPLEMENT_CONST_MOCK0(getName);
};

#endif //MBL_ED_ANALYZERTASKMOCK_H
