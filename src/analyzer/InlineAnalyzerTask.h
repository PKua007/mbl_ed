//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_INLINEANALYZERTASK_H
#define MBL_ED_INLINEANALYZERTASK_H


#include "AnalyzerTask.h"

/**
 * @brief An AnalyzerTask, whose result is in the form of a few numbers (or strings) with labels.
 *
 * Multiple InlineAnalyzerTask -s, together with additional simulation parameters are supposed to form a single row
 * in a results table, see InlineResultsPrinter.
 */
class InlineAnalyzerTask : public AnalyzerTask {
public:
    /**
     * @brief Returns a vector with descriptions of result fields with the same order as in getResultFields().
     */
    [[nodiscard]] virtual std::vector<std::string> getResultHeader() const = 0;

    /**
     * @brief Returns a vector with values of result fields with the same order as in getResultHeader().
     */
    [[nodiscard]] virtual std::vector<std::string> getResultFields() const = 0;
};


#endif //MBL_ED_INLINEANALYZERTASK_H
