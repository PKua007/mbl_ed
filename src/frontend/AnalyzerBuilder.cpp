//
// Created by pkua on 08.06.2020.
//

#include <memory>

#include "AnalyzerBuilder.h"

#include "utils/Assertions.h"
#include "analyzer/tasks/CDF.h"
#include "analyzer/tasks/MeanGapRatio.h"
#include "analyzer/tasks/MeanInverseParticipationRatio.h"
#include "analyzer/tasks/InverseParticipationRatio.h"
#include "analyzer/tasks/EDCorrelationsTimeEvolution.h"

Analyzer AnalyzerBuilder::build(const std::vector<std::string> &tasks, const Parameters &params,
                                                 std::shared_ptr<FockBase> fockBase) {
    Analyzer analyzer;
    for (const auto &task : tasks) {
        std::istringstream taskStream(task);
        std::string taskName;
        taskStream >> taskName;
        if (taskName == "mgr") {
            double mgrCenter, mgrMargin;
            taskStream >> mgrCenter >> mgrMargin;
            ValidateMsg(taskStream, "Wrong format, use: mgr [epsilon center] [epsilon margin]");
            Validate(mgrCenter > 0 && mgrCenter < 1);
            Validate(mgrMargin > 0 && mgrMargin <= 1);
            Validate(mgrCenter - mgrMargin/2 >= 0 && mgrCenter + mgrMargin/2 <= 1);
            analyzer.addTask(std::make_unique<MeanGapRatio>(mgrCenter, mgrMargin));
        } else if (taskName == "mipr") {
            double mgrCenter, mgrMargin;
            taskStream >> mgrCenter >> mgrMargin;
            ValidateMsg(taskStream, "Wrong format, use: mipr [epsilon center] [epsilon margin]");
            Validate(mgrCenter > 0 && mgrCenter < 1);
            Validate(mgrMargin > 0 && mgrMargin <= 1);
            Validate(mgrCenter - mgrMargin/2 >= 0 && mgrCenter + mgrMargin/2 <= 1);
            analyzer.addTask(std::make_unique<MeanInverseParticipationRatio>(mgrCenter, mgrMargin));
        } else if (taskName == "ipr") {
            double mgrCenter, mgrMargin;
            taskStream >> mgrCenter >> mgrMargin;
            ValidateMsg(taskStream, "Wrong format, use: ipr [epsilon center] [epsilon margin]");
            Validate(mgrCenter > 0 && mgrCenter < 1);
            Validate(mgrMargin > 0 && mgrMargin <= 1);
            Validate(mgrCenter - mgrMargin/2 >= 0 && mgrCenter + mgrMargin/2 <= 1);
            analyzer.addTask(std::make_unique<InverseParticipationRatio>(mgrCenter, mgrMargin));
        } else if (taskName == "cdf") {
            std::size_t bins;
            taskStream >> bins;
            ValidateMsg(taskStream, "Wrong format, use: cdf [number of bins]");
            Validate(bins >= 2);
            analyzer.addTask(std::make_unique<CDF>(bins));
        } else if (taskName == "evolution") {
            if (params.N != params.K || params.K % 2 != 0)
                throw ValidationException("evolution mode is only for even number of sites with 1:1 filling");

            CorrelationsTimeEvolutionParameters evolutionParameters;
            evolutionParameters.fockBase = fockBase;
            evolutionParameters.numberOfSites = params.K;

            std::string vectorsToEvolveStr;
            double maxTime;
            std::size_t numSteps;
            taskStream >> maxTime >> numSteps >> evolutionParameters.marginSize;
            evolutionParameters.timeSegmentation.push_back({maxTime, numSteps});
            taskStream >> vectorsToEvolveStr;
            ValidateMsg(taskStream, "Wrong format, use: evolution [max time] [number of steps] [margin size] "
                                    "[vectors to evolve - unif/dw/both]\nunif - 1.1.1.1; dw - 2.0.2.0; both - both ;)");
            Validate(evolutionParameters.timeSegmentation[0].maxTime > 0);
            Validate(evolutionParameters.timeSegmentation[0].numSteps >= 2);
            Validate(evolutionParameters.marginSize * 2 < params.K);

            evolutionParameters.setVectorsToEvolveFromTag(vectorsToEvolveStr);

            analyzer.addTask(std::make_unique<EDCorrelationsTimeEvolution>(evolutionParameters));
        } else {
            throw ValidationException("Unknown analyzer task: " + taskName);
        }
    }
    return analyzer;
}
