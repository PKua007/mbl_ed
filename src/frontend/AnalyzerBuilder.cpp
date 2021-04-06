//
// Created by pkua on 08.06.2020.
//

#include <memory>
#include <iterator>
#include <algorithm>
#include <ZipIterator.hpp>

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
    // All parsing facilities are global functions, not (static) methods, because we have AnalyzerBuilder mainly for
    // a "compilation wall" in the Frontend - it doesn't have to be compiled after every alteration

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

    BandExtractor::Margin parse_margin(std::istringstream &taskStream, std::size_t basisDim, const std::string &usage) {
        BandExtractor::Margin margin;

        std::string secondArg{};
        taskStream >> secondArg;
        ValidateMsg(taskStream, usage);

        if (secondArg == "n") {
            std::size_t numEnergies{};
            taskStream >> numEnergies;
            ValidateMsg(taskStream, usage);
            Validate(numEnergies > 0);
            Validate(numEnergies <= basisDim);
            margin = BandExtractor::NumberOfEnergiesMargin(numEnergies);
        } else {
            double epsilonMargin{parse_double(secondArg)};
            Validate(epsilonMargin > 0 && epsilonMargin <= 1);
            margin = BandExtractor::WidthMargin(epsilonMargin);
        }
        return margin;
    }

    FockVector parse_fock_vector(std::size_t numSites, const std::string &vectorString, const std::string &usage) {
        FockBasis::Vector middleVector;

        try {
            // Try form dw|unif [epsilon margin]
            middleVector = FockBasis::Vector(numSites, vectorString);
        } catch (FockVectorParseException &) {
            // Last chance - 1.2.0.0 [epsilon margin], otherwise error
            try {
                middleVector = FockBasis::Vector(vectorString);
            } catch (std::invalid_argument &) {
                throw ValidationException(usage);
            }
        }
        return middleVector;
    }

    BandExtractor::Range parse_band(std::size_t numSites, std::size_t basisDim, std::istringstream &taskStream,
                                    const std::string &taskName)
    {
        BandExtractor::Range range = BandExtractor::EpsilonRange(0.5, 0.1);

        std::string firstArg;
        taskStream >> firstArg;
        std::string usage = "Wrong format, use one of:\n"
                            + taskName + " [epsilon center] [epsilon margin]|{n [number of energies]}\n"    // case 1
                            + taskName + " dw|unif|2.1.0.0 [epsilon margin]|{n [number of energies]}\n"     // case 2
                            + taskName + " cdf [cdf center] [cdf margin]";                                  // case 3
        ValidateMsg(taskStream, usage);

        if (firstArg == "cdf") {
            // case 3
            double cdfMiddle{}, cdfMargin{};
            taskStream >> cdfMiddle >> cdfMargin;
            ValidateMsg(taskStream, usage);
            Validate(cdfMiddle > 0 && cdfMiddle < 1);
            Validate(cdfMargin > 0 && cdfMargin <= 1);
            Validate(cdfMiddle - cdfMargin / 2 >= 0 && cdfMiddle + cdfMargin / 2 <= 1);
            range = BandExtractor::CDFRange(cdfMiddle, cdfMargin);
        } else {
            // case 1 or case 2
            BandExtractor::Margin margin = parse_margin(taskStream, basisDim, usage);

            const std::string &middleEnergyMarker = firstArg;
            try {
                // First, try form [epsilon center] [epsilon margin]
                double epsilonMiddle = parse_double(middleEnergyMarker);
                Validate(epsilonMiddle > 0 && epsilonMiddle < 1);
                range = BandExtractor::EpsilonRange(epsilonMiddle, margin);
            } catch (ParseDoubleException &) {
                // If not, try parsing vector
                FockVector middleVector = parse_fock_vector(numSites, middleEnergyMarker, usage);
                range = BandExtractor::VectorRange(middleVector, margin);
            }
        }
        return range;
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
            auto band = parse_band(params.K, fockBasis->size(), taskStream, "mgr");
            analyzer->addTask(std::make_unique<MeanGapRatio>(band));
        } else if (taskName == "mipr") {
            auto band = parse_band(params.K, fockBasis->size(), taskStream, "mipr");
            analyzer->addTask(std::make_unique<MeanInverseParticipationRatio>(band));
        } else if (taskName == "ipr") {
            auto band = parse_band(params.K, fockBasis->size(), taskStream, "ipr");
            analyzer->addTask(std::make_unique<InverseParticipationRatio>(band));
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

            std::array<std::string, 3> fields = {std::string(" time "), std::string(" obs "), std::string(" vec ")};
            std::array<std::size_t, 3> idxs{};

            for (std::size_t i{}; i < fields.size(); i++) {
                const auto &field = fields[i];
                std::size_t idx = task.find(field);
                ValidateMsg(idx != std::string::npos, "Malformed evolution params - missing field: " + field);
                idxs[i] = idx;
            }

            auto zipped = Zip(idxs, fields);
            std::sort(zipped.begin(), zipped.end());

            std::array<std::string, 3> fieldValues{};
            for (std::size_t i{}; i < fields.size(); i++) {
                const auto &field = fields[i];
                auto idx = idxs[i];
                auto &fieldValue = fieldValues[i];

                std::size_t beg = idx + field.length();
                if (i == fields.size() - 1) {
                    ValidateMsg(beg < task.length(), "Missing field value for: " + field);
                    fieldValue = task.substr(beg);
                    trim(fieldValue);
                } else {
                    std::size_t end = idxs[i + 1];
                    ValidateMsg(beg < end, "Missing gield value for: " + field);
                    fieldValue = task.substr(beg, end - beg);
                    trim(fieldValue);
                }
                Validate(!fieldValue.empty());
            }

            TimeEvolutionParameters evolutionParams;
            evolutionParams.fockBasis = fockBasis;
            evolutionParams.numberOfSites = params.K;
            ObservablesBuilder observablesBuilder;

            for (std::size_t i{}; i < fields.size(); i++) {
                const auto &field = fields[i];
                const auto &fieldValue = fieldValues[i];

                if (field == " time ") {
                    std::istringstream timeSegmentationStream(fieldValue);
                    std::copy(std::istream_iterator<EvolutionTimeSegment>(timeSegmentationStream),
                              std::istream_iterator<EvolutionTimeSegment>(),
                              std::back_inserter(evolutionParams.timeSegmentation));
                } else if (field == " obs ") {
                    std::vector<std::string> observablesParams = explode(fieldValue, ';');
                    observablesBuilder.build(observablesParams, params, fockBasis, hamiltonianGenerator);
                } else if (field == " vec ") {
                    std::istringstream tagsStream(fieldValue);
                    std::vector<std::string> tags;
                    std::copy(std::istream_iterator<std::string>(tagsStream), std::istream_iterator<std::string>(),
                              std::back_inserter(tags));
                    evolutionParams.setVectorsToEvolveFromTags(tags);
                } else {
                    throw std::runtime_error(field);
                }
            }

            evolutionParams.primaryObservables = observablesBuilder.releasePrimaryObservables();
            evolutionParams.secondaryObservables = observablesBuilder.releaseSecondaryObservables();
            evolutionParams.storedObservables = observablesBuilder.releaseStoredObservables();

            auto observablesEvolution = std::make_unique<OservablesTimeEvolution>();
            analyzer->addTask(std::make_unique<EDTimeEvolution>(evolutionParams, std::move(observablesEvolution)));
        } else if (taskName == "dressed") {
            double coeffThreshold;
            taskStream >> coeffThreshold;
            Validate(coeffThreshold > M_SQRT1_2);
            auto band = parse_band(params.K, fockBasis->size(), taskStream,
                                   "dressed [coefficient threshold > 1/sqrt(2)]");
            analyzer->addTask(std::make_unique<DressedStatesFinder>(coeffThreshold, band));
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
            double q;
            taskStream >> q;
            Validate(q > 0);
            auto band = parse_band(params.K, fockBasis->size(), taskStream,"pe [q parameter]");
            analyzer->addTask(std::make_unique<ParticipationEntropy>(q, band));
        } else {
            throw ValidationException("Unknown analyzer task: " + taskName);
        }
    }
    return analyzer;
}
