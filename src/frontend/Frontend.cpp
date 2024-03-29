//
// Created by Piotr Kubala on 21/01/2020.
//

#include <cxxopts.hpp>
#include <filesystem>

#include "Frontend.h"
#include "HamiltonianGeneratorBuilder.h"
#include "AnalyzerBuilder.h"
#include "AveragingModelFactory.h"
#include "IO.h"
#include "ObservablesBuilder.h"

#include "simulation/ExactDiagonalization.h"
#include "core/FockBasisGenerator.h"
#include "core/DisorderGenerator.h"

#include "core/QuenchCalculator.h"
#include "simulation/ChebyshevEvolution.h"
#include "evolution/TimeEvolution.h"

#include "simulation/RestorableSimulationExecutor.h"
#include "simulation/QuenchDataSimulation.h"

#include "utils/Fold.h"
#include "utils/Utils.h"
#include "utils/Assertions.h"
#include "utils/OMPMacros.h"
#include "simulation/RandomStateObservables.h"


void Frontend::ed(int argc, char **argv) {
    // Parse options
    cxxopts::Options options(argv[0],
                             Fold("Exact diagonalization mode. Performs one or more simulations depending on the "
                             "parameters. It can store them for to further processing and/or perform some "
                             "analyzer tasks on the fly.")
                             .width(80));

    std::string inputFilename;
    std::filesystem::path outputFilename;
    std::filesystem::path directory;
    std::vector<std::string> paramsToPrint;
    std::vector<std::string> onTheFlyTasks;
    std::vector<std::string> overridenParams;
    std::string verbosity;

    options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters. See input.ini for parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("o,output", "when specified, inline results will be printed to this file",
             cxxopts::value<std::filesystem::path>(outputFilename))
            ("t,task", "task(s) to be performed on the fly while simulating",
             cxxopts::value<std::vector<std::string>>(onTheFlyTasks))
            ("d,directory", "where to put eigenenergies and auxiliary analyzer files (like per-eigensystem observable"
                            "values). Note, that this option does not affect --output path nor bulk analyzer result "
                            "files",
             cxxopts::value<std::filesystem::path>(directory)->default_value("."))
            ("p,print_parameter", "parameters to be included in inline results",
             cxxopts::value<std::vector<std::string>>(paramsToPrint)
                 ->default_value("N,K"))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P N=1 (-PN=1 does not work) act as one would append N=1 to [general] section of input"
                            "file. To override or even add some hamiltonian terms use -P termName.paramName=value",
             cxxopts::value<std::vector<std::string>>(overridenParams))
            ("V,verbosity", "how verbose the output should be. Allowed values, with increasing verbosity: "
                            "error, warn, info, verbose, debug",
             cxxopts::value<std::string>(verbosity)->default_value("info"));

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        this->out << options.help() << std::endl;
        exit(0);
    }

    Logger logger(this->out);
    this->setOverridenParamsAsAdditionalText(logger, overridenParams);
    this->setVerbosityLevel(logger, verbosity);

    // Validate parsed options
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help", logger);
    if (!parsedOptions.count("input"))
        die("Input file must be specified with option -i [input file name]", logger);
    if (!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory))
        die("Output directory " + directory.string() + " does not exist or is not a directory", logger);

    // Prepare parameters
    IO io(logger);
    Parameters params = io.loadParameters(inputFilename, overridenParams);
    params.print(logger);
    logger << std::endl;
    for (const auto &param : paramsToPrint)
        if (!params.hasParam(param))
            die("Parameters to print: parameter " + param + " is unknown", logger);

    // OpenMP info
    logger.info() << "Using " << _OMP_MAXTHREADS << " OpenMP threads" << std::endl;

    // Generate Fock basis
    FockBasisGenerator basisGenerator;
    logger.verbose() << "Preparing Fock basis started... " << std::endl;
    arma::wall_clock timer;
    timer.tic();
    auto basis = std::shared_ptr<FockBasis>(basisGenerator.generate(params.N, params.K));
    logger.info() << "Preparing Fock basis done (" << timer.toc() << " s)." << std::endl;

    // Prepare HamiltonianGenerator, Analyzer and AveragingModel
    auto rnd = std::make_unique<RND>();
    auto hamiltonianGenerator = HamiltonianGeneratorBuilder{}.build(params, basis, *rnd);
    auto analyzer = AnalyzerBuilder{}.build(onTheFlyTasks, params, basis, *hamiltonianGenerator, directory);
    auto averagingModel = AveragingModelFactory{}.create(params.averagingModel);

    // Prepare and run simulations
    ExactDiagonalizationParameters simulationParams = this->prepareExactDiagonalizationParameters(directory, params);
    ExactDiagonalization simulation(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                    simulationParams, std::move(analyzer));

    SimulationsSpan simulationsSpan;
    simulationsSpan.from = params.from;
    simulationsSpan.to = params.to;
    simulationsSpan.total = params.totalSimulations;
    RestorableSimulationExecutor restorableSimulationExecutor(simulationsSpan, params.getOutputFileSignatureWithRange(),
                                                              params.splitWorkload, params.secureSimulationState);
    restorableSimulationExecutor.performSimulations(simulation, params.seed, logger);

    // Save results
    if (restorableSimulationExecutor.shouldSaveSimulation()) {
        const Analyzer &analyzerRef = simulation.getAnalyzer();
        io.printInlineResults(params, paramsToPrint, analyzerRef.getInlineResultsHeader(),
                              analyzerRef.getInlineResultsFields());
        if (outputFilename.empty())
            io.storeAnalyzerResults(params, analyzerRef, paramsToPrint, std::nullopt);
        else
            io.storeAnalyzerResults(params, analyzerRef, paramsToPrint, outputFilename);
    }
}

ExactDiagonalizationParameters Frontend::prepareExactDiagonalizationParameters(const std::filesystem::path &directory,
                                                                               const Parameters &params) const
{
    ExactDiagonalizationParameters simulationParams;
    simulationParams.calculateEigenvectors = params.calculateEigenvectors;

    if (!params.saveEigenenergies && !params.saveEigenstates)
        simulationParams.storeLevel = ExactDiagonalizationParameters::StoreLevel::NONE;
    else if (params.saveEigenenergies && !params.saveEigenstates)
        simulationParams.storeLevel = ExactDiagonalizationParameters::StoreLevel::EIGENENERGIES;
    else if (params.saveEigenenergies && params.saveEigenstates)
        simulationParams.storeLevel = ExactDiagonalizationParameters::StoreLevel::EIGENSYSTEM;
    else
        throw std::runtime_error("");

    if (params.storeFormat.empty() || params.storeFormat == "arma_binary" || params.storeFormat == "default")
        simulationParams.fileType = arma::arma_binary;
    else if (params.storeFormat == "arma_ascii")
        simulationParams.fileType = arma::arma_ascii;
    else if (params.storeFormat == "raw_ascii")
        simulationParams.fileType = arma::raw_ascii;
    else if (params.storeFormat == "csv_ascii")
        simulationParams.fileType = arma::csv_ascii;
    else
        throw ValidationException("Unknown store format: " + params.storeFormat);

    simulationParams.fileSignature = directory / params.getOutputFileSignature();
    return simulationParams;
}

void Frontend::setOverridenParamsAsAdditionalText(Logger &logger, std::vector<std::string> overridenParams) const {
    if (overridenParams.empty())
        return;

    for (auto &param : overridenParams) {
        std::regex pattern(R"([a-zA-Z0-9_]+\.(.+=.+))");
        std::string replaced = std::regex_replace(param, pattern, "$1");
        if (!replaced.empty())
            param = replaced;
    }

    std::ostringstream additionalText;
    std::copy(overridenParams.begin(), overridenParams.end() - 1,
              std::ostream_iterator<std::string>(additionalText, ", "));
    additionalText << overridenParams.back();
    logger.setAdditionalText(additionalText.str());
}

void Frontend::analyze(int argc, char **argv) {
    // Parse options
    cxxopts::Options options(argv[0],
                             Fold("Mode for analyzing the results of already performed simulations, stored in the "
                             "files. It performs one or more analyzer task, specified with option -t (--task). The "
                             "files are found using signature generated from parameters.").width(80));

    std::string inputFilename;
    std::string outputFilename;
    std::vector<std::string> paramsToPrint{};
    std::vector<std::string> tasks{};
    std::filesystem::path directory;
    std::vector<std::string> overridenParams;
    std::string verbosity;

    options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters. See input.ini for parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("o,output", "when specified, inline results will be printed to this file",
             cxxopts::value<std::string>(outputFilename))
            ("t,task", "analyzer task(s) to be performed",
             cxxopts::value<std::vector<std::string>>(tasks))
            ("p,print_parameter", "parameters to be included in inline results",
             cxxopts::value<std::vector<std::string>>(paramsToPrint)->default_value("N,K"))
            ("d,directory", "directory to search simulation results",
             cxxopts::value<std::filesystem::path>(directory)->default_value("."))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P N=1 (-PN=1 does not work) act as one would append N=1 to [general] section of input"
                            "file. To override or even add some hamiltonian terms use -P termName.paramName=value",
             cxxopts::value<std::vector<std::string>>(overridenParams))
            ("V,verbosity", "how verbose the output should be. Allowed values, with increasing verbosity: "
                            "error, warn, info, verbose, debug",
             cxxopts::value<std::string>(verbosity)->default_value("info"));

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        this->out << options.help() << std::endl;
        exit(0);
    }

    Logger logger(this->out);
    this->setOverridenParamsAsAdditionalText(logger, overridenParams);
    this->setVerbosityLevel(logger, verbosity);

    // Validate options
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help", logger);
    if (!parsedOptions.count("input"))
        die("Input file must be specified with option -i [input file name]", logger);
    if (!parsedOptions.count("task"))
        die("At least 1 analyzer task must be specified with option -t [task parameters]", logger);
    if (!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory))
        die("Directory " + directory.string() + " does not exist or is not a directory", logger);

    // Load parameters
    IO io(logger);
    Parameters params = io.loadParameters(inputFilename, overridenParams);
    params.print(logger);
    logger << std::endl;
    for (const auto &param : paramsToPrint)
        if (!params.hasParam(param))
            die("Parameters to print: parameter " + param + " is unknown", logger);

    // Generate Fock basis
    FockBasisGenerator basisGenerator;
    logger.verbose() << "Preparing Fock basis started... " << std::endl;
    arma::wall_clock timer;
    timer.tic();
    auto basis = std::shared_ptr<FockBasis>(basisGenerator.generate(params.N, params.K));
    logger.info() << "Preparing Fock basis done (" << timer.toc() << " s)." << std::endl;

    // Load eigenenergies and analyze them
    auto analyzer = AnalyzerBuilder{}.build(tasks, params, basis, std::nullopt, directory);
    std::string fileSignature = params.getOutputFileSignature();
    std::vector<std::string> energiesFilenames = io.findFiles(directory, fileSignature, "nrg.bin");
    std::vector<std::string> statesFilenames = io.findFiles(directory, fileSignature, "st.bin");
    this->validateEigensystemFiles(logger, energiesFilenames, statesFilenames);

    for (const auto &energiesFilename : energiesFilenames) {
        bool restoreEigenstates = !statesFilenames.empty();
        Eigensystem eigensystem = this->restoreEigensystem(energiesFilename, restoreEigenstates, basis, logger);
        logger.verbose() << "Analyzing " << energiesFilename << " started... " << std::endl;
        timer.tic();
        analyzer->analyze(eigensystem, logger);
        logger.info() << "Analyzing " << energiesFilename << " done (" << timer.toc() << " s)." << std::endl;
    }

    // Save results
    io.printInlineResults(params, paramsToPrint, analyzer->getInlineResultsHeader(), analyzer->getInlineResultsFields());
    if (outputFilename.empty())
        io.storeAnalyzerResults(params, *analyzer, paramsToPrint, std::nullopt);
    else
        io.storeAnalyzerResults(params, *analyzer, paramsToPrint, outputFilename);
}

Eigensystem Frontend::restoreEigensystem(const std::string &energiesFilename, bool restoreEigenstates,
                                         std::shared_ptr<FockBasis> basis, Logger &logger) const
{
    std::ifstream energiesFile(energiesFilename);
    if (!energiesFile)
        die("Cannot open " + energiesFilename + " to read eigenenergies from", logger);

    Eigensystem eigensystem;
    if (restoreEigenstates) {
        auto statesFilename = this->stripSuffix(energiesFilename, "nrg.bin") + "st.bin";
        std::ifstream statesFile(statesFilename);
        if (!energiesFile)
            die("Cannot open " + statesFilename + " to read eigenstates from", logger);

        eigensystem.restore(energiesFile, statesFile, basis);
    } else {
        eigensystem.restore(energiesFile, basis);
    }

    return eigensystem;
}

void Frontend::validateEigensystemFiles(Logger &logger, const std::vector<std::string> &energiesFilenames,
                                        const std::vector<std::string> &statesFilenames)
{
    if (energiesFilenames.empty()) {
        die("No eigenenergy files were found.", logger);
    } else if (!statesFilenames.empty()) {
        std::vector<std::string> notMatching;
        this->appendNotMatchingSignatures(notMatching, energiesFilenames, "nrg.bin", statesFilenames, "st.bin");
        this->appendNotMatchingSignatures(notMatching, statesFilenames, "st.bin", energiesFilenames, "nrg.bin");

        if (!notMatching.empty()) {
            logger.error() << "Both eigenenergy and eigenstate files were found, but the signatures of the following ";
            logger << "files don't match:" << std::endl;
            for (const auto &signature : notMatching)
                logger << ">> " << signature << std::endl;
            die("");
        }

        logger.info() << "Found " << energiesFilenames.size() << " eigenenergy and eigenstate files" << std::endl;
    } else {
        logger.info() << "Found " << energiesFilenames.size() << " eigenenergy files" << std::endl;
    }
}

void Frontend::chebyshev(int argc, char **argv) {
    // Parse options
    cxxopts::Options options(argv[0],
                             Fold("Performs evolution using Chebyshev expansion technique.").width(80));

    std::string inputFilename;
    std::vector<std::string> overridenParamsEntries;
    std::vector<std::string> quenchParamsEntries;
    std::string timeSegmentation{};
    std::vector<std::string> vectorsToEvolveTags;
    std::string verbosity;
    std::vector<std::string> observableStrings;

    options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters. See input.ini for parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P N=1 (-PN=1 does not work) act as one would append N=1 to [general] section of input "
                            "file. To override or even add some hamiltonian terms use -P termName.paramName=value",
             cxxopts::value<std::vector<std::string>>(overridenParamsEntries))
            ("t,time_segmentation", "describes the time span and how it should be divided, in format: [time 1] [number "
                                    "of steps 1] [time 2] [number of steps 2] ... . For example, '1 2 5 4' divides 0-1 "
                                    "in 2 and 1-5 in 4, giving times: 0, 0.5, 1, 2, 3, 4, 5",
             cxxopts::value<std::string>(timeSegmentation))
            ("O,observable", "which observables should be calculated and stored. One or more of: "
                             "n_j - onsite occupations; "
                             "n_iN_j - products of 2 onsite occupations; "
                             "G_d [margin size] - site averaged correlations, without including [margin size] sites on "
                             "both ends; "
                             "rho_i - fluctiations;",
             cxxopts::value<std::vector<std::string>>(observableStrings)->default_value("G_d 0, rho_i, n_i"))
            ("v,vectors", "vectors to evolve. Available options: unif/dw/2.3.0.0.1. You can specify more than one",
             cxxopts::value<std::vector<std::string>>(vectorsToEvolveTags))
            ("q,quench_param", "if specified, evolution will be performed for quenched vector (together with "
                                "those specified by --vectors). This option overrides the param as in --set_param, "
                                "applied after --set_param, but for the separate initial Hamiltonian in quantum quench",
             cxxopts::value<std::vector<std::string>>(quenchParamsEntries))
            ("V,verbosity", "how verbose the output should be. Allowed values, with increasing verbosity: "
                            "error, warn, info, verbose, debug",
             cxxopts::value<std::string>(verbosity)->default_value("info"));

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        this->out << options.help() << std::endl;
        exit(0);
    }

    Logger logger(this->out);
    std::vector allOverridenEntries(overridenParamsEntries);
    allOverridenEntries.insert(allOverridenEntries.end(), quenchParamsEntries.begin(), quenchParamsEntries.end());
    this->setOverridenParamsAsAdditionalText(logger, allOverridenEntries);
    this->setVerbosityLevel(logger, verbosity);

    // Validate parsed options - input file
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help", logger);
    if (!parsedOptions.count("input"))
        die("Input file must be specified with option -i [input file name]", logger);

    // Prepare parameters
    IO io(logger);
    Parameters params = io.loadParameters(inputFilename, overridenParamsEntries);
    params.print(logger);
    logger << std::endl;

    // Validate rest of the options
    if (!parsedOptions.count("time_segmentation"))
        die("You have to specify max evolution time using -t [max time]", logger);
    if (!parsedOptions.count("vectors") && !parsedOptions.count("quench_param"))
        die("You have to specify space vectors to evolve using -v [unif/dw/1.0.4.0] or/and via quench -q", logger);
    // Validation of vectors is done later

    // Prepare quench parameters, if desired
    std::optional<Parameters> quenchParams;
    if (parsedOptions.count("quench_param")) {
        std::vector<std::string> overridenAndQuenchParamEntries = overridenParamsEntries;
        overridenAndQuenchParamEntries.insert(overridenAndQuenchParamEntries.end(), quenchParamsEntries.begin(),
                                              quenchParamsEntries.end());
        quenchParams = io.loadParameters(inputFilename, overridenAndQuenchParamEntries);

        logger.info() << "Together with vectors specified by -v, evolution will be performed for a state quenched ";
        logger << "from initial Hamiltonian:" << std::endl;
        quenchParams->printHamiltonianTerms(logger);
        logger << std::endl;
    }

    // OpenMP info
    logger.info() << "Using " << _OMP_MAXTHREADS << " OpenMP threads" << std::endl;

    // Prepare FockBasis
    FockBasisGenerator basisGenerator;
    logger.verbose() << "Preparing Fock basis started... " << std::endl;
    arma::wall_clock timer;
    timer.tic();
    auto basis = std::shared_ptr<FockBasis>(basisGenerator.generate(params.N, params.K));
    logger.info() << "Preparing Fock basis done (" << timer.toc() << " s)." << std::endl;

    // Prepare HamiltonianGenerator, Analyzer and AveragingModel
    auto rnd = std::make_unique<RND>();
    auto hamiltonianGenerator = HamiltonianGeneratorBuilder{}.build(params, basis, *rnd);
    auto averagingModel = AveragingModelFactory{}.create(params.averagingModel);

    // Prepare and run evolutions
    TimeEvolutionParameters evolutionParams;

    std::istringstream timeSegmentationStream(timeSegmentation);
    std::copy(std::istream_iterator<EvolutionTimeSegment>(timeSegmentationStream),
              std::istream_iterator<EvolutionTimeSegment>(),
              std::back_inserter(evolutionParams.timeSegmentation));
    evolutionParams.numberOfSites = params.K;
    evolutionParams.fockBasis = basis;
    evolutionParams.setVectorsToEvolveFromTags(vectorsToEvolveTags); // This one also does the validation

    ObservablesBuilder observablesBuilder;
    observablesBuilder.build(observableStrings, params, basis, *hamiltonianGenerator);
    evolutionParams.primaryObservables = observablesBuilder.releasePrimaryObservables();
    evolutionParams.secondaryObservables = observablesBuilder.releaseSecondaryObservables();
    evolutionParams.storedObservables = observablesBuilder.releaseStoredObservables();

    auto observablesEvolution = std::make_unique<OservablesTimeEvolution>();

    SimulationsSpan simulationsSpan;
    simulationsSpan.from = params.from;
    simulationsSpan.to = params.to;
    simulationsSpan.total = params.totalSimulations;
    RestorableSimulationExecutor simulationExecutor(simulationsSpan, params.getOutputFileSignatureWithRange(),
                                                    params.splitWorkload, params.secureSimulationState);

    std::unique_ptr<ChebyshevEvolution<>> evolution;
    if (quenchParams.has_value()) {
        using ExternalVector = TimeEvolutionParameters::ExternalVector;
        evolutionParams.vectorsToEvolve.emplace_back(ExternalVector{"quench"});

        auto quenchRnd = std::make_unique<RND>();
        auto quenchHamiltonianGenerator = HamiltonianGeneratorBuilder{}.build(*quenchParams, basis, *quenchRnd);
        evolution = std::make_unique<ChebyshevEvolution<>>(
                std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                std::make_unique<TimeEvolution>(evolutionParams, std::move(observablesEvolution)),
                std::make_unique<QuenchCalculator>(), std::move(quenchHamiltonianGenerator), std::move(quenchRnd)
        );
    } else {
        evolution = std::make_unique<ChebyshevEvolution<>>(
            std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
            std::make_unique<TimeEvolution>(evolutionParams, std::move(observablesEvolution))
        );
    }

    simulationExecutor.performSimulations(*evolution, params.seed, logger);
    evolution->printQuenchInfo(logger);

    if (simulationExecutor.shouldSaveSimulation()) {
        std::string resultsFilename = params.getOutputFileSignatureWithRange() + "_evolution.txt";
        std::ofstream resultsFile(resultsFilename);
        evolution->storeResults(resultsFile);
        logger.info() << "Observables stored to " << resultsFilename << std::endl;
    }
}

void Frontend::quench(int argc, char **argv) {
    // Parse options
    cxxopts::Options options(argv[0], Fold("Performs quantum quench by fining the ground state of initial "
                                           "Hamiltonian and calculating its energy and quantum spread in final "
                                           "Hamiltonian, where some parameters are different. The values are averaged "
                                           "over multiple realisations.").width(80));

    std::string inputFilename;
    std::string outputFilename;
    std::vector<std::string> overridenParamsEntries;
    std::vector<std::string> quenchParamsEntries;
    std::vector<std::string> paramsToPrint{};
    std::string verbosity;

    options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters. See input.ini for parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("o,output", "when specified, quench results will be printed to this file",
             cxxopts::value<std::string>(outputFilename))
            ("p,print_parameter", "parameters to be included in inline results",
             cxxopts::value<std::vector<std::string>>(paramsToPrint)->default_value("N,K"))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P N=1 (-PN=1 does not work) act as one would append N=1 to [general] section of input"
                            "file. To override or even add some hamiltonian terms use -P termName.paramName=value",
             cxxopts::value<std::vector<std::string>>(overridenParamsEntries))
            ("q,quench_param", "overrides the param as in --set_param, applied after --set_param, but for the initial"
                               "Hamiltonian in quantum quench",
             cxxopts::value<std::vector<std::string>>(quenchParamsEntries))
            ("V,verbosity", "how verbose the output should be. Allowed values, with increasing verbosity: "
                            "error, warn, info, verbose, debug",
             cxxopts::value<std::string>(verbosity)->default_value("info"));

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        this->out << options.help() << std::endl;
        exit(0);
    }

    Logger logger(this->out);
    std::vector allOverridenEntries(overridenParamsEntries);
    allOverridenEntries.insert(allOverridenEntries.end(), quenchParamsEntries.begin(), quenchParamsEntries.end());
    this->setOverridenParamsAsAdditionalText(logger, allOverridenEntries);
    this->setVerbosityLevel(logger, verbosity);

    // Validate parsed options
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help", logger);
    if (!parsedOptions.count("input"))
        die("Input file must be specified with option -i [input file name]", logger);
    if (quenchParamsEntries.empty())
        die("At least one parameter should be specified for quantum quench", logger);

    // Prepare quench (initial) parameters
    IO io(logger);
    std::vector<std::string> overridenAndQuenchParamEntries = overridenParamsEntries;
    overridenAndQuenchParamEntries.insert(overridenAndQuenchParamEntries.end(), quenchParamsEntries.begin(),
                                          quenchParamsEntries.end());
    Parameters quenchParams = io.loadParameters(inputFilename, overridenAndQuenchParamEntries);
    for (const auto &param : paramsToPrint)
        if (!quenchParams.hasParam(param))
            die("Parameters to print: parameter " + param + " is unknown", logger);

    // Prepare (final) parameters
    Parameters params = io.loadParameters(inputFilename, overridenParamsEntries);

    // Print info on prameters and initial and final Hamiltonians
    params.printGeneral(logger);
    logger << std::endl;
    logger.info() << "Initial Hamiltonian:" << std::endl;
    quenchParams.printHamiltonianTerms(logger);
    logger << std::endl;
    logger.info() << "Final Hamiltonian:" << std::endl;
    params.printHamiltonianTerms(logger);
    logger << std::endl;

    // Generate Fock basis
    FockBasisGenerator basisGenerator;
    logger.verbose() << "Preparing Fock basis started... " << std::endl;
    arma::wall_clock timer;
    timer.tic();
    auto basis = std::shared_ptr<FockBasis>(basisGenerator.generate(params.N, params.K));
    logger.info() << "Preparing Fock basis done (" << timer.toc() << " s)." << std::endl;

    // Prepare initial and final HamiltonianGenerator
    auto initialRnd = std::make_unique<RND>();
    auto finalRnd = std::make_unique<RND>();
    auto initialHamiltonianGenerator = HamiltonianGeneratorBuilder{}.build(quenchParams, basis, *initialRnd);
    auto finalHamiltonianGenerator = HamiltonianGeneratorBuilder{}.build(params, basis, *finalRnd);
    auto averagingModel = AveragingModelFactory{}.create(params.averagingModel);

    SimulationsSpan simulationsSpan;
    simulationsSpan.from = params.from;
    simulationsSpan.to = params.to;
    simulationsSpan.total = params.totalSimulations;
    RestorableSimulationExecutor simulationExecutor(simulationsSpan, params.getOutputFileSignatureWithRange(),
                                                    params.splitWorkload, params.secureSimulationState);

    // Prepare and run quenches
    auto quenchCalculator = std::make_unique<QuenchCalculator>();
    QuenchDataSimulation simulation(std::move(initialHamiltonianGenerator), std::move(initialRnd),
                                    std::move(finalHamiltonianGenerator), std::move(finalRnd),
                                    std::move(averagingModel), std::move(quenchCalculator));
    simulationExecutor.performSimulations(simulation, params.seed, logger);

    // Save results
    if (simulationExecutor.shouldSaveSimulation()) {
        std::vector<std::string> resultHeader = simulation.getResultsHeader();
        std::vector<std::string> resultFields = simulation.getResultsFields();
        io.printInlineResults(quenchParams, paramsToPrint, resultHeader, resultFields);
        if (!outputFilename.empty())
            io.storeInlineResults(quenchParams, paramsToPrint, resultHeader, resultFields, outputFilename);
    }
}

void Frontend::randomStates(int argc, char **argv) {
    // Parse options
    cxxopts::Options options(argv[0], "Mode calculating mean observable values for random Gaussian states.");

    std::string inputFilename;
    std::filesystem::path outputFilename;
    std::vector<std::string> paramsToPrint;
    std::vector<std::string> overridenParams;
    std::vector<std::string> observableStrings;
    std::string verbosity;

    options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters. See input.ini for parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("o,output", "when specified, results will be printed to this file",
             cxxopts::value<std::filesystem::path>(outputFilename))
            ("O,observable", "which observables should be calculated and stored. One or more of: "
                             "n_j - onsite occupations; "
                             "n_iN_j - products of 2 onsite occupations; "
                             "G_d [margin size] - site averaged correlations, without including [margin size] sites on "
                             "both ends; "
                             "rho_i - fluctiations;",
             cxxopts::value<std::vector<std::string>>(observableStrings)->default_value("G_d 0, rho_i, n_i"))
            ("p,print_parameter", "parameters to be included in inline results",
             cxxopts::value<std::vector<std::string>>(paramsToPrint)->default_value("N,K"))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P N=1 (-PN=1 does not work) act as one would append N=1 to [general] section of input"
                            "file. To override or even add some hamiltonian terms use -P termName.paramName=value",
             cxxopts::value<std::vector<std::string>>(overridenParams))
            ("V,verbosity", "how verbose the output should be. Allowed values, with increasing verbosity: "
                            "error, warn, info, verbose, debug",
             cxxopts::value<std::string>(verbosity)->default_value("info"));

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        this->out << options.help() << std::endl;
        exit(0);
    }

    Logger logger(this->out);
    this->setOverridenParamsAsAdditionalText(logger, overridenParams);
    this->setVerbosityLevel(logger, verbosity);

    // Validate parsed options
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help", logger);
    if (!parsedOptions.count("input"))
        die("Input file must be specified with option -i [input file name]", logger);

    // Prepare parameters
    IO io(logger);
    Parameters params = io.loadParameters(inputFilename, overridenParams);
    params.print(logger);
    logger << std::endl;
    for (const auto &param : paramsToPrint)
        if (!params.hasParam(param))
            die("Parameters to print: parameter " + param + " is unknown", logger);

    // OpenMP info
    logger.info() << "Using " << _OMP_MAXTHREADS << " OpenMP threads" << std::endl;

    // Generate Fock basis
    FockBasisGenerator basisGenerator;
    logger.verbose() << "Preparing Fock basis started... " << std::endl;
    arma::wall_clock timer;
    timer.tic();
    auto basis = std::shared_ptr<FockBasis>(basisGenerator.generate(params.N, params.K));
    logger.info() << "Preparing Fock basis done (" << timer.toc() << " s)." << std::endl;

    // Prepare simulation classes
    auto rnd = std::make_unique<RND>();

    ObservablesBuilder observablesBuilder;
    observablesBuilder.build(observableStrings, params, basis, std::nullopt);
    RandomStateObservables simulation(std::move(rnd), std::move(basis), observablesBuilder.releasePrimaryObservables(),
                                      observablesBuilder.releaseSecondaryObservables(),
                                      observablesBuilder.releaseStoredObservables());

    SimulationsSpan simulationsSpan;
    simulationsSpan.from = params.from;
    simulationsSpan.to = params.to;
    simulationsSpan.total = params.totalSimulations;
    RestorableSimulationExecutor restorableSimulationExecutor(simulationsSpan, params.getOutputFileSignatureWithRange(),
                                                              params.splitWorkload, params.secureSimulationState);
    restorableSimulationExecutor.performSimulations(simulation, params.seed, logger);

    // Save results
    if (restorableSimulationExecutor.shouldSaveSimulation()) {
        io.printInlineResults(params, paramsToPrint, simulation.getHeader(), simulation.getValues());
        if (!outputFilename.empty())
            io.storeInlineResults(params, paramsToPrint, simulation.getHeader(), simulation.getValues(), outputFilename);
    }
}

void Frontend::preview(int argc, char **argv) {
    // Parse options
    cxxopts::Options options(argv[0], "Auxiliary mode to preview different stuff.");

    std::string inputFilename;
    std::vector<std::string> overridenParams;
    std::string verbosity;

    options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters, needed when some other options are specified. See input.ini for "
                        "parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("P,set-param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P N=1 (-PN=1 does not work) act as one would append N=1 to [general] section of the "
                            "input file. To override or even add some hamiltonian terms use -P termName.paramName=value",
             cxxopts::value<std::vector<std::string>>(overridenParams))
            ("F,preview-fock-basis", "when specified, Fock basis will be printed (in raw format). --input option needs "
                                     "to be specified")
            ("V,verbosity", "how verbose the output should be. Allowed values, with increasing verbosity: "
                            "error, warn, info, verbose, debug",
             cxxopts::value<std::string>(verbosity)->default_value("info"))
            ("c,clear-output", "only the proper output is printed, without logging info, except for errors. "
                               "Overrides --verbosity.");

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        this->out << options.help() << std::endl;
        exit(0);
    }

    Logger logger(this->out);
    this->setOverridenParamsAsAdditionalText(logger, overridenParams);
    this->setVerbosityLevel(logger, verbosity);
    if (parsedOptions.count("clear-output"))
        logger.setVerbosityLevel(Logger::LogType::ERROR);

    // Validate parsed options
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help", logger);

    // Load parameters (if specified)
    std::optional<Parameters> params;
    if (parsedOptions.count("input")) {
        IO io(logger);
        params = io.loadParameters(inputFilename, overridenParams);
        if (!parsedOptions.count("clear-output")) {
            params->print(logger);
            logger << std::endl;
        }
    }

    // Fock basis preview mode
    if (parsedOptions.count("preview-fock-basis")) {
        if (!parsedOptions.count("input"))
            die("To preview the Fock basis, --input option must be specified");
        Assert(params.has_value());

        FockBasisGenerator basisGenerator;
        logger.verbose() << "Preparing Fock basis started... " << std::endl;
        arma::wall_clock timer;
        timer.tic();
        auto basis = std::shared_ptr<FockBasis>(basisGenerator.generate(params->N, params->K));
        logger.info() << "Preparing Fock basis done (" << timer.toc() << " s)." << std::endl;

        std::ostream &rawLogger = logger;
        for (const auto &fockVector : *basis)
            rawLogger << fockVector << std::endl;

        return;
    }

    die("No preview mode specified. See " + cmd + " --help");
}

void Frontend::printGeneralHelp(const std::string &cmd) {
    this->out << Fold("Program performing exact diagonalization of Hubbard-like hamiltonians with analyzing "
                      "facilities. ").width(80) << std::endl;
    this->out << std::endl;
    this->out << "Usage: " << cmd << " [mode] (mode dependent parameters). " << std::endl;
    this->out << std::endl;
    this->out << "Available modes:" << std::endl;
    this->out << "ed" << std::endl;
    this->out << Fold("Performs the diagonalizations according to the parameters passed. Some analyzer tasks can "
                      "be performed on the fly.").width(80).margin(4) << std::endl;
    this->out << "analyze" << std::endl;
    this->out << Fold("Performs one or more analyzer tasks after loading simulation results from the files.")
                 .width(80).margin(4) << std::endl;
    this->out << "chebyshev" << std::endl;
    this->out << Fold("Performs time evolution using Chebyshev expansion technique.")
                 .width(80).margin(4) << std::endl;
    this->out << "quench" << std::endl;
    this->out << Fold("Performs quantum quench from some initial to final Hamiltonian and print energy info.")
                 .width(80).margin(4) << std::endl;
    this->out << "random-states" << std::endl;
    this->out << Fold("Calculates mean observable values for random Gaussian states.")
                 .width(80).margin(4) << std::endl;
    this->out << "preview" << std::endl;
    this->out << Fold("Auxiliary mode to preview different stuff.")
                 .width(80).margin(4) << std::endl;
    this->out << std::endl;
    this->out << "Type " + cmd + " [mode] --help to get help on the specific mode." << std::endl;
}

void Frontend::setVerbosityLevel(Logger &logger, const std::string &verbosityLevelName) const {
    if (verbosityLevelName == "error")
        logger.setVerbosityLevel(Logger::ERROR);
    else if (verbosityLevelName == "warn")
        logger.setVerbosityLevel(Logger::WARN);
    else if (verbosityLevelName == "info")
        logger.setVerbosityLevel(Logger::INFO);
    else if (verbosityLevelName == "verbose")
        logger.setVerbosityLevel(Logger::VERBOSE);
    else if (verbosityLevelName == "debug")
        logger.setVerbosityLevel(Logger::DEBUG);
    else
        die("Unknown verbosity level: " + verbosityLevelName, logger);
}

void Frontend::appendNotMatchingSignatures(std::vector<std::string> &notMatchingSignatures,
                                           const std::vector<std::string> &filenames1, const std::string &suffix1,
                                           const std::vector<std::string> &filenames2, const std::string &suffix2)
{
    for (const auto &filename1 : filenames1) {
        Assert(endsWith(filename1, suffix1));
        std::string filenameSuffix2 = this->stripSuffix(filename1, suffix1) + suffix2;
        if (std::find(filenames2.begin(), filenames2.end(), filenameSuffix2) == filenames2.end())
            notMatchingSignatures.push_back(filename1);
    }
}

std::string Frontend::stripSuffix(const std::string &string, const std::string &suffix) const {
    Assert(endsWith(string, suffix));
    return string.substr(0, string.length() - suffix.length());
}
