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
#include "analyzer/tasks/BulkMeanGapRatio.h"
#include "evolution/observables/OnsiteOccupations.h"
#include "evolution/observables/OnsiteOccupationsSquared.h"
#include "evolution/observables/Correlations.h"
#include "evolution/observables/OnsiteFluctuations.h"


std::unique_ptr<Analyzer> AnalyzerBuilder::build(const std::vector<std::string> &tasks, const Parameters &params,
                                                 std::shared_ptr<FockBase> fockBase)
{
    auto analyzer = std::make_unique<Analyzer>();
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
                analyzer->addTask(std::make_unique<MeanGapRatio>(*middleVector, mgrMargin));
            } else {
                double mgrCenter = std::stod(mgrCenterString);
                Validate(mgrCenter > 0 && mgrCenter < 1);
                Validate(mgrMargin > 0 && mgrMargin <= 1);
                Validate(mgrCenter - mgrMargin / 2 >= 0 && mgrCenter + mgrMargin / 2 <= 1);
                analyzer->addTask(std::make_unique<MeanGapRatio>(mgrCenter, mgrMargin));
            }
        } else if (taskName == "mipr") {
            double mgrCenter, mgrMargin;
            taskStream >> mgrCenter >> mgrMargin;
            ValidateMsg(taskStream, "Wrong format, use: mipr [epsilon center] [epsilon margin]");
            Validate(mgrCenter > 0 && mgrCenter < 1);
            Validate(mgrMargin > 0 && mgrMargin <= 1);
            Validate(mgrCenter - mgrMargin/2 >= 0 && mgrCenter + mgrMargin/2 <= 1);
            analyzer->addTask(std::make_unique<MeanInverseParticipationRatio>(mgrCenter, mgrMargin));
        } else if (taskName == "ipr") {
            double mgrCenter, mgrMargin;
            taskStream >> mgrCenter >> mgrMargin;
            ValidateMsg(taskStream, "Wrong format, use: ipr [epsilon center] [epsilon margin]");
            Validate(mgrCenter > 0 && mgrCenter < 1);
            Validate(mgrMargin > 0 && mgrMargin <= 1);
            Validate(mgrCenter - mgrMargin/2 >= 0 && mgrCenter + mgrMargin/2 <= 1);
            analyzer->addTask(std::make_unique<InverseParticipationRatio>(mgrCenter, mgrMargin));
        } else if (taskName == "cdf") {
            std::size_t bins;
            taskStream >> bins;
            ValidateMsg(taskStream, "Wrong format, use: cdf [number of bins]");
            Validate(bins >= 2);
            analyzer->addTask(std::make_unique<CDF>(bins));
        } else if (taskName == "evolution") {
            if (params.N != params.K || params.K % 2 != 0)
                throw ValidationException("evolution mode is only for even number of sites with 1:1 filling");

            TimeEvolutionParameters evolutionParams;
            evolutionParams.fockBase = fockBase;
            evolutionParams.numberOfSites = params.K;

            double maxTime;
            std::size_t numSteps;
            std::size_t marginSize;
            taskStream >> maxTime >> numSteps >> marginSize;
            evolutionParams.timeSegmentation.emplace_back(EvolutionTimeSegment{maxTime, numSteps});
            ValidateMsg(taskStream, "Wrong format, use: evolution [max time] [number of steps] [margin size] "
                                    "[space-separated vectors to evolve - unif/dw/1.0.4.0]");
            std::vector<std::string> tags;
            std::copy(std::istream_iterator<std::string>(taskStream), std::istream_iterator<std::string>(),
                      std::back_inserter(tags));
            Validate(evolutionParams.timeSegmentation[0].maxTime > 0);
            Validate(evolutionParams.timeSegmentation[0].numSteps >= 2);
            Validate(marginSize * 2 < params.K);

            evolutionParams.setVectorsToEvolveFromTags(tags);

            auto occupations = std::make_shared<OnsiteOccupations>(fockBase);
            auto occupationsSquared = std::make_shared<OnsiteOccupationsSquared>(fockBase);
            auto correlations = std::make_shared<Correlations>(params.K, 0);
            auto borderlessCorrelations = std::make_shared<Correlations>(params.K, marginSize);
            auto fluctuations = std::make_shared<OnsiteFluctuations>(params.K);
            evolutionParams.primaryObservables = {occupations, occupationsSquared};
            evolutionParams.secondaryObservables = {correlations, borderlessCorrelations, fluctuations};
            evolutionParams.storedObservables = {correlations, borderlessCorrelations, fluctuations, occupations};

            auto occupationEvolution = std::make_unique<OservablesTimeEvolution>(
                evolutionParams.primaryObservables, evolutionParams.secondaryObservables,
                evolutionParams.storedObservables
            );

            analyzer->addTask(
                std::make_unique<EDCorrelationsTimeEvolution>(evolutionParams, std::move(occupationEvolution))
            );
        } else if (taskName == "dressed") {
            double mgrCenter, mgrMargin, coeffThreshold;
            taskStream >> mgrCenter >> mgrMargin >> coeffThreshold;
            ValidateMsg(taskStream, "Wrong format, use: dressed [epsilon center] [epsilon margin] "
                                    "[coefficient threshold > 1/sqrt(2)");
            Validate(mgrCenter > 0 && mgrCenter < 1);
            Validate(mgrMargin > 0 && mgrMargin <= 1);
            Validate(mgrCenter - mgrMargin / 2 >= 0 && mgrCenter + mgrMargin / 2 <= 1);
            Validate(coeffThreshold > M_SQRT1_2);
            analyzer->addTask(std::make_unique<DressedStatesFinder>(mgrCenter, mgrMargin, coeffThreshold));
        } else if (taskName == "mgrs") {
            std::size_t numBins;
            taskStream >> numBins;
            ValidateMsg(taskStream, "Wrong format, use: mgrs [number of bins]");
            Validate(numBins > 0);
            analyzer->addTask(std::make_unique<BulkMeanGapRatio>(numBins));
        } else {
            throw ValidationException("Unknown analyzer task: " + taskName);
        }
    }
    return analyzer;
}
