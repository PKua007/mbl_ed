//
// Created by pkua on 08.06.2020.
//

#include <memory>
#include <iterator>
#include <algorithm>

#include "AnalyzerBuilder.h"

#include "utils/Assertions.h"
#include "analyzer/tasks/CDF.h"
#include "analyzer/tasks/MeanGapRatio.h"
#include "analyzer/tasks/MeanInverseParticipationRatio.h"
#include "analyzer/tasks/InverseParticipationRatio.h"
#include "analyzer/tasks/EDCorrelationsTimeEvolution.h"
#include "analyzer/tasks/DressedStatesFinder.h"

Analyzer AnalyzerBuilder::build(const std::vector<std::string> &tasks, const Parameters &params,
                                std::shared_ptr<FockBase> fockBase)
{
    Analyzer analyzer;
    for (const auto &task : tasks) {
        std::istringstream taskStream(task);
        std::string taskName;
        taskStream >> taskName;
        if (taskName == "mgr") {
            std::string mgrCenterString;
            double mgrMargin;
            taskStream >> mgrCenterString >> mgrMargin;
            ValidateMsg(taskStream, "Wrong format, use: mgr [epsilon center - number or 'unif' or 'dw']"
                                    " [epsilon margin]");

            std::optional<FockBase::Vector> middleVector;
            if (mgrCenterString == "unif") {
                ValidateMsg(params.K == params.N, "'unif' and 'dw' are only for a unity filling");
                middleVector = FockBase::Vector(params.K, 1);
            } else if (mgrCenterString == "dw") {
                ValidateMsg(params.K == params.N, "'unif' and 'dw' are only for a unity filling");
                ValidateMsg(params.K % 2 == 0, "'dw' is only available for even number of sites");
                middleVector = FockBase::Vector(params.K);
                for (std::size_t i{}; i < middleVector->size(); i += 2)
                    (*middleVector)[i] = 2;
            }

            if (middleVector.has_value()) {
                analyzer.addTask(std::make_unique<MeanGapRatio>(*middleVector, mgrMargin));
            } else {
                double mgrCenter = std::stod(mgrCenterString);
                Validate(mgrCenter > 0 && mgrCenter < 1);
                Validate(mgrMargin > 0 && mgrMargin <= 1);
                Validate(mgrCenter - mgrMargin / 2 >= 0 && mgrCenter + mgrMargin / 2 <= 1);
                analyzer.addTask(std::make_unique<MeanGapRatio>(mgrCenter, mgrMargin));
            }
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

            double maxTime;
            std::size_t numSteps;
            taskStream >> maxTime >> numSteps >> evolutionParameters.marginSize;
            evolutionParameters.timeSegmentation.push_back({maxTime, numSteps});
            ValidateMsg(taskStream, "Wrong format, use: evolution [max time] [number of steps] [margin size] "
                                    "[space-separated vectors to evolve - unif/dw/1.0.4.0]");
            std::vector<std::string> tags;
            std::copy(std::istream_iterator<std::string>(taskStream), std::istream_iterator<std::string>(),
                      std::back_inserter(tags));
            Validate(evolutionParameters.timeSegmentation[0].maxTime > 0);
            Validate(evolutionParameters.timeSegmentation[0].numSteps >= 2);
            Validate(evolutionParameters.marginSize * 2 < params.K);

            evolutionParameters.setVectorsToEvolveFromTag(tags);

            analyzer.addTask(std::make_unique<EDCorrelationsTimeEvolution>(evolutionParameters));
        } else if (taskName == "dressed") {
            double mgrCenter, mgrMargin, coeffThreshold;
            taskStream >> mgrCenter >> mgrMargin >> coeffThreshold;
            ValidateMsg(taskStream, "Wrong format, use: dressed [epsilon center] [epsilon margin] "
                                    "[coefficient threshold > 1/sqrt(2)");
            Validate(mgrCenter > 0 && mgrCenter < 1);
            Validate(mgrMargin > 0 && mgrMargin <= 1);
            Validate(mgrCenter - mgrMargin / 2 >= 0 && mgrCenter + mgrMargin / 2 <= 1);
            Validate(coeffThreshold > M_SQRT1_2);
            analyzer.addTask(std::make_unique<DressedStatesFinder>(mgrCenter, mgrMargin, coeffThreshold));
        } else {
            throw ValidationException("Unknown analyzer task: " + taskName);
        }
    }
    return analyzer;
}
