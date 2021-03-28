//
// Created by pkua on 08.06.2020.
//

#include <memory>
#include <iterator>
#include <algorithm>

#include "AnalyzerBuilder.h"
#include "ObservablesBuilder.h"

#include "utils/Assertions.h"
#include "utils/Utils.h"

#include "analyzer/tasks/CDF.h"
#include "analyzer/tasks/PDF.h"
#include "analyzer/tasks/MeanGapRatio.h"
#include "analyzer/tasks/MeanInverseParticipationRatio.h"
#include "analyzer/tasks/InverseParticipationRatio.h"
#include "analyzer/tasks/EDTimeEvolution.h"
#include "analyzer/tasks/DressedStatesFinder.h"
#include "analyzer/tasks/BulkMeanGapRatio.h"
#include "analyzer/tasks/EigenstateObservables.h"
#include "analyzer/tasks/ParticipationEntropy.h"

#include "core/observables/OnsiteOccupations.h"
#include "core/observables/OnsiteOccupationsSquared.h"
#include "core/observables/Correlations.h"
#include "core/observables/OnsiteFluctuations.h"


namespace {
    struct ParseDoubleException : public std::runtime_error {
        explicit ParseDoubleException(const std::string &what) : std::runtime_error(what) { }
    };

    double parse_double(std::string str) {
        try {
            trim(str);
            std::size_t offset = 0;
            double result = std::stod(str, &offset);
            if (offset != str.size())
                throw ParseDoubleException("Additional characters found while parsing double");
            return result;
        } catch (const std::invalid_argument &e) {
            throw ParseDoubleException(e.what());
        }
    }
}

std::unique_ptr<Analyzer> AnalyzerBuilder::build(const std::vector<std::string> &tasks, const Parameters &params,
                                                 const std::shared_ptr<FockBasis> &fockBasis,
                                                 std::optional<std::reference_wrapper<const HamiltonianGenerator>>
                                                 hamiltonianGenerator, const std::filesystem::path &auxiliaryDir)
{
    auto analyzer = std::make_unique<Analyzer>();
    for (const auto &task : tasks) {
        std::istringstream taskStream(task);
        std::string taskName;
        taskStream >> taskName;
        if (taskName == "mgr") {
            std::string firstArg;
            taskStream >> firstArg;
            std::string usage = "Wrong format, use one of:\n"
                                "mgr [epsilon center] [epsilon margin]|{n [number of energies]}\n"
                                "mgr dw|unif|2.1.0.0 [epsilon margin]|{n [number of energies]}\n"
                                "mgr cdf [cdf center] [cdf margin]";
            ValidateMsg(taskStream, usage);

            if (firstArg == "cdf") {
                double cdfMiddle{}, cdfMargin{};
                taskStream >> cdfMiddle >> cdfMargin;
                ValidateMsg(taskStream, usage);
                Validate(cdfMiddle > 0 && cdfMiddle < 1);
                Validate(cdfMargin > 0 && cdfMargin <= 1);
                Validate(cdfMiddle - cdfMargin / 2 >= 0 && cdfMiddle + cdfMargin / 2 <= 1);
                analyzer->addTask(std::make_unique<MeanGapRatio>(BandExtractor::CDFRange(cdfMiddle, cdfMargin)));
            } else {
                const std::string &middle = firstArg;
                std::string marginToken{};
                taskStream >> marginToken;
                ValidateMsg(taskStream, usage);

                double epsilonMargin{};
                std::size_t numEnergies{};

                if (marginToken == "n") {
                    taskStream >> numEnergies;
                    ValidateMsg(taskStream, usage);
                    Validate(numEnergies > 0);
                    Validate(numEnergies <= fockBasis->size());
                } else {
                    epsilonMargin = parse_double(marginToken);
                    Validate(epsilonMargin > 0 && epsilonMargin <= 1);
                }

                // First, try form [epsilon center] [epsilon margin]
                try {
                    double epsilonMiddle = parse_double(middle);
                    Validate(epsilonMiddle > 0 && epsilonMiddle < 1);

                    if (marginToken == "n") {
                        analyzer->addTask(
                            std::make_unique<MeanGapRatio>(BandExtractor::EpsilonRange(epsilonMiddle, numEnergies))
                        );
                    } else {
                        Validate(epsilonMiddle - epsilonMargin / 2 >= 0 && epsilonMiddle + epsilonMargin / 2 <= 1);
                        analyzer->addTask(
                            std::make_unique<MeanGapRatio>(BandExtractor::EpsilonRange(epsilonMiddle, epsilonMargin))
                        );
                    }
                } catch (ParseDoubleException &) {
                    FockBasis::Vector middleVector;

                    // Next, try form dw|unif [epsilon margin]
                    try {
                        middleVector = FockBasis::Vector(params.K, middle);
                    } catch (FockVectorParseException &) {
                        // Last chance - 1.2.0.0 [epsilon margin], otherwise error
                        try {
                            middleVector = FockBasis::Vector(middle);
                        } catch (std::invalid_argument &) {
                            throw ValidationException(usage);
                        }
                    }

                    if (marginToken == "n") {
                        analyzer->addTask(
                            std::make_unique<MeanGapRatio>(BandExtractor::VectorRange(middleVector, numEnergies))
                        );
                    } else {
                        analyzer->addTask(
                            std::make_unique<MeanGapRatio>(BandExtractor::VectorRange(middleVector, epsilonMargin))
                        );
                    }
                }
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
        } else if (taskName == "pdf") {
            std::size_t bins;
            taskStream >> bins;
            ValidateMsg(taskStream, "Wrong format, use: pdf [number of bins]");
            Validate(bins >= 2);
            analyzer->addTask(std::make_unique<PDF>(bins));
        } else if (taskName == "evolution") {
            if (params.N != params.K || params.K % 2 != 0)
                throw ValidationException("evolution mode is only for even number of sites with 1:1 filling");

            TimeEvolutionParameters evolutionParams;
            evolutionParams.fockBasis = fockBasis;
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

            auto occupations = std::make_shared<OnsiteOccupations>(fockBasis);
            auto occupationsSquared = std::make_shared<OnsiteOccupationsSquared>(fockBasis);
            auto correlations = std::make_shared<Correlations>(params.K, 0);
            auto borderlessCorrelations = std::make_shared<Correlations>(params.K, marginSize);
            auto fluctuations = std::make_shared<OnsiteFluctuations>(params.K);
            evolutionParams.primaryObservables = {occupations, occupationsSquared};
            evolutionParams.secondaryObservables = {correlations, borderlessCorrelations, fluctuations};
            evolutionParams.storedObservables = {correlations, borderlessCorrelations, fluctuations, occupations};

            auto occupationEvolution = std::make_unique<OservablesTimeEvolution>();

            analyzer->addTask(
                std::make_unique<EDTimeEvolution>(evolutionParams, std::move(occupationEvolution))
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
        } else if (taskName == "obs") {
            bool store{};
            if (taskStream.str().find("store") != std::string::npos) {
                std::string storeString;
                taskStream >> storeString;
                ValidateMsg(taskStream && storeString == "store", "Wrong format, use: obs (store) [number of bins] "
                                                                  "[obs. 1] [obs. 1 params] [obs. 2] "
                                                                  "[obs. 2 params];...");
                store = true;
            }
            std::size_t numBins;
            taskStream >> numBins;
            ValidateMsg(taskStream, "Wrong format, use: obs (store) [number of bins] [obs. 1] [obs. 1 params];"
                                    "[obs. 2] [obs. 2 params];...");
            std::string observablesString;
            std::getline(taskStream, observablesString);
            std::vector<std::string> observablesParams = explode(observablesString, ';');
            ObservablesBuilder builder;
            builder.build(observablesParams, params, fockBasis, hamiltonianGenerator);
            auto eigenstateObservables = std::make_unique<EigenstateObservables>(
                numBins, builder.releasePrimaryObservables(), builder.releaseSecondaryObservables(),
                builder.releaseStoredObservables()
            );
            if (store)
                eigenstateObservables->startStoringObservables(auxiliaryDir / params.getOutputFileSignature());
            analyzer->addTask(std::move(eigenstateObservables));
        } else if (taskName == "pe") {
            double q, mgrCenter, mgrMargin;
            taskStream >> q >> mgrCenter >> mgrMargin;
            ValidateMsg(taskStream, "Wrong format, use: pe [q parameter] [epsilon center] [epsilon margin]");
            Validate(mgrCenter > 0 && mgrCenter < 1);
            Validate(mgrMargin > 0 && mgrMargin <= 1);
            Validate(mgrCenter - mgrMargin/2 >= 0 && mgrCenter + mgrMargin/2 <= 1);
            Validate(q > 0);
            analyzer->addTask(std::make_unique<ParticipationEntropy>(q, mgrCenter, mgrMargin));
        } else {
            throw ValidationException("Unknown analyzer task: " + taskName);
        }
    }
    return analyzer;
}
