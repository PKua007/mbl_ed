//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_INLINEANALYZERTASKMOCK_H
#define MBL_ED_INLINEANALYZERTASKMOCK_H

#include <catch2/trompeloeil.hpp>

#include "analyzer/InlineAnalyzerTask.h"

class InlineAnalyzerTaskMock  : public trompeloeil::mock_interface<InlineAnalyzerTask> {
public:
    IMPLEMENT_CONST_MOCK0(getResultHeader);
    IMPLEMENT_CONST_MOCK0(getResultFields);
    IMPLEMENT_MOCK1(analyze);
    IMPLEMENT_CONST_MOCK0(getName);
};


#endif //MBL_ED_INLINEANALYZERTASKMOCK_H
