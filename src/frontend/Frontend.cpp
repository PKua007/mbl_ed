//
// Created by Piotr Kubala on 21/01/2020.
//

#include <cxxopts.hpp>
#include <filesystem>
#include <omp.h>

#include "Frontend.h"
#include "HamiltonianGeneratorBuilder.h"
#include "AnalyzerBuilder.h"
#include "AveragingModelFactory.h"
#include "IO.h"

#include "simulation/Simulation.h"
#include "simulation/FockBaseGenerator.h"
#include "simulation/DisorderGenerator.h"
#include "simulation/QuenchCalculator.h"

#include "evolution/ChebyshevEvolution.h"

#include "utils/Fold.h"
#include "utils/Utils.h"
#include "utils/Assertions.h"


void Frontend::simulate(int argc, char **argv) {
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

    options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters. See input.ini for parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("o,output", "when specified, inline results will be printed to this file",
             cxxopts::value<std::filesystem::path>(outputFilename))
            ("t,task", "task(s) to be performed on the fly while simulating",
             cxxopts::value<std::vector<std::string>>(onTheFlyTasks))
            ("d,directory", "where to put eigenenergies and eigenstates. Note, that this option does not "
                            "affect --output path nor analyzer result files",
             cxxopts::value<std::filesystem::path>(directory)->default_value("."))
            ("p,print_parameter", "parameters to be included in inline results",
             cxxopts::value<std::vector<std::string>>(paramsToPrint)
                 ->default_value("N,K"))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P N=1 (-PN=1 does not work) act as one would append N=1 to [general] section of input"
                            "file. To override or even add some hamiltonian terms use -P termName.paramName=value",
             cxxopts::value<std::vector<std::string>>(overridenParams));

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Validate parsed options
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help");
    if (!parsedOptions.count("input"))
        die("Input file must be specified with option -i [input file name]");
    if (!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory))
        die("Output directory " + directory.string() + " does not exist or is not a directory");

    // Prepare parameters
    IO io(std::cout);
    Parameters params = io.loadParameters(inputFilename, overridenParams);
    params.print(std::cout);
    std::cout << std::endl;
    for (const auto &param : paramsToPrint)
        if (!params.hasParam(param))
            die("Parameters to print: parameter " + param + " is unknown");

    // OpenMP info
    std::cout << "[Frontend::simulate] Using " << omp_get_max_threads() << " OpenMP threads" << std::endl;

    // Generate Fock basis
    FockBaseGenerator baseGenerator;
    std::cout << "[Frontend::simulate] Preparing Fock basis... " << std::flush;
    arma::wall_clock timer;
    timer.tic();
    auto base = std::shared_ptr(baseGenerator.generate(params.N, params.K));
    std::cout << "done (" << timer.toc() << " s)." << std::endl;

    // Prepare HamiltonianGenerator, Analyzer and AveragingModel
    auto rnd = std::make_unique<RND>(params.from + params.seed);
    auto hamiltonianGenerator = HamiltonianGeneratorBuilder{}.build(params, base, *rnd);
    Analyzer analyzer = AnalyzerBuilder{}.build(onTheFlyTasks, params, base);
    auto averagingModel = AveragingModelFactory{}.create(params.averagingModel);

    // Prepare and run simulations
    SimulationParameters simulationParams;
    simulationParams.from = params.from;
    simulationParams.to = params.to;
    simulationParams.totalSimulations = params.totalSimulations;
    simulationParams.calculateEigenvectors = params.calculateEigenvectors;
    simulationParams.saveEigenenergies = params.saveEigenenergies;
    simulationParams.fileSignature = directory / params.getOutputFileSignature();

    Simulation simulation(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd), simulationParams);
    simulation.perform(this->out, analyzer);

    // Save results
    io.printInlineResults(params, paramsToPrint, analyzer.getInlineResultsHeader(), analyzer.getInlineResultsFields());
    if (outputFilename.empty())
        io.storeAnalyzerResults(params, analyzer, paramsToPrint, std::nullopt);
    else
        io.storeAnalyzerResults(params, analyzer, paramsToPrint, outputFilename);
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

    options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters. See input.ini for parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("o,output", "when specified, inline results will be printed to this file",
             cxxopts::value<std::string>(outputFilename))
            ("t,task", "analyzer task(s) to be performed",
             cxxopts::value<std::vector<std::string>>(tasks))
            ("p,print_parameter", "parameters to be included in inline results",
             cxxopts::value<std::vector<std::string>>(paramsToPrint)
                 ->default_value("N,K"))
            ("d,directory", "directory to search simulation results",
             cxxopts::value<std::filesystem::path>(directory)->default_value("."))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P N=1 (-PN=1 does not work) act as one would append N=1 to [general] section of input"
                            "file. To override or even add some hamiltonian terms use -P termName.paramName=value",
             cxxopts::value<std::vector<std::string>>(overridenParams));

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Validate options
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help");
    if (!parsedOptions.count("input"))
        die("Input file must be specified with option -i [input file name]");
    if (!parsedOptions.count("task"))
        die("At least 1 analyzer task must be specified with option -t [task parameters]");
    if (!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory))
        die("Directory " + directory.string() + " does not exist or is not a directory");

    // Load parameters
    IO io(std::cout);
    Parameters params = io.loadParameters(inputFilename, overridenParams);
    params.print(std::cout);
    std::cout << std::endl;
    for (const auto &param : paramsToPrint)
        if (!params.hasParam(param))
            die("Parameters to print: parameter " + param + " is unknown");

    // Generate Fock basis
    FockBaseGenerator baseGenerator;
    std::cout << "[Frontend::analyze] Preparing Fock basis... " << std::flush;
    arma::wall_clock timer;
    timer.tic();
    auto base = std::shared_ptr(baseGenerator.generate(params.N, params.K));
    std::cout << "done (" << timer.toc() << " s)." << std::endl;

    // Load eigenenergies and analyze them
    Analyzer analyzer = AnalyzerBuilder{}.build(tasks, params, base);
    std::string fileSignature = params.getOutputFileSignature();
    std::vector<std::string> energiesFilenames = io.findEigenenergyFiles(directory, fileSignature);
    if (energiesFilenames.empty())
        die("No eigenenergy files were found.");

    for (const auto &energiesFilename : energiesFilenames) {
        std::ifstream energiesFile(energiesFilename);
        if (!energiesFile)
            die("Cannot open " + energiesFilename + " to read eigenenergies from");

        Eigensystem eigensystem;
        eigensystem.restore(energiesFile, base);
        analyzer.analyze(eigensystem, std::cout);
    }

    // Save results
    io.printInlineResults(params, paramsToPrint, analyzer.getInlineResultsHeader(), analyzer.getInlineResultsFields());
    if (outputFilename.empty())
        io.storeAnalyzerResults(params, analyzer, paramsToPrint, std::nullopt);
    else
        io.storeAnalyzerResults(params, analyzer, paramsToPrint, outputFilename);
}

void Frontend::chebyshev(int argc, char **argv) {
    // Parse options
    cxxopts::Options options(argv[0],
                             Fold("Performs evolution using Chebyshev expansion technique.").width(80));

    std::string inputFilename;
    std::vector<std::string> overridenParams;
    std::string timeSegmentation{};
    std::size_t marginSize{};
    std::vector<std::string> vectorsToEvolveTags;

    options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters. See input.ini for parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P N=1 (-PN=1 does not work) act as one would append N=1 to [general] section of input "
                            "file. To override or even add some hamiltonian terms use -P termName.paramName=value",
             cxxopts::value<std::vector<std::string>>(overridenParams))
            ("t,time_segmentation", "describes the time span and how it should be divided, in format: [time 1] [number "
                                    "of steps 1] [time 2] [number of steps 2] ... . For example, '1 2 5 4' divides 0-1 "
                                    "in 2 and 1-5 in 4, giving times: 0, 0.5, 1, 2, 3, 4, 5",
             cxxopts::value<std::string>(timeSegmentation))
            ("m,margin", "margin size - one averaging of correlations is done for all sites, the second one for all "
                         "but a given margin from both sides", cxxopts::value<std::size_t>(marginSize))
            ("v,vectors", "vectors to evolve. Available options: unif/dw/both; unif - 1.1.1.1; dw - 2.0.2.0;"
                          " both - both ;)", cxxopts::value<std::vector<std::string>>(vectorsToEvolveTags));

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Validate parsed options - input file
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help");
    if (!parsedOptions.count("input"))
        die("Input file must be specified with option -i [input file name]");

    // Prepare parameters
    IO io(std::cout);
    Parameters params = io.loadParameters(inputFilename, overridenParams);
    params.print(std::cout);
    std::cout << std::endl;

    // Validate rest of the options
    if (!parsedOptions.count("time_segmentation"))
        die("You have to specify max evolution time using -t [max time]");
    if (!parsedOptions.count("margin"))
        die("You have to specify margin size using -m [margin size]");
    if (marginSize * 2 > params.K - 2)
        die("Margin is too big - there should be at least 2 sites left.");
    if (!parsedOptions.count("vectors"))
        die("You have to specify space vectors to evolve using -v [unif/dw/1.0.4.0]");
    // Validation of vectors is done later

    // OpenMP info
    std::cout << "[Frontend::chebyshev] Using " << omp_get_max_threads() << " OpenMP threads" << std::endl;

    // Prepare FockBasis
    FockBaseGenerator baseGenerator;
    std::cout << "[Frontend::chebyshev] Preparing Fock basis... " << std::flush;
    arma::wall_clock timer;
    timer.tic();
    auto base = std::shared_ptr(baseGenerator.generate(params.N, params.K));
    std::cout << "done (" << timer.toc() << " s)." << std::endl;

    // Prepare HamiltonianGenerator, Analyzer and AveragingModel
    auto rnd = std::make_unique<RND>(params.from + params.seed);
    auto hamiltonianGenerator = HamiltonianGeneratorBuilder{}.build(params, base, *rnd);
    auto averagingModel = AveragingModelFactory{}.create(params.averagingModel);

    // Prepare and run evolutions
    CorrelationsTimeEvolutionParameters evolutionParameters;

    std::istringstream timeSegmentationStream(timeSegmentation);
    std::copy(std::istream_iterator<EvolutionTimeSegment>(timeSegmentationStream),
              std::istream_iterator<EvolutionTimeSegment>(),
              std::back_inserter(evolutionParameters.timeSegmentation));
    evolutionParameters.numberOfSites = params.K;
    evolutionParameters.fockBase = base;
    evolutionParameters.marginSize = marginSize;
    evolutionParameters.setVectorsToEvolveFromTags(vectorsToEvolveTags); // This one also does the validation

    ChebyshevEvolution evolution(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                 params.from, params.to, params.totalSimulations, evolutionParameters,
                                 params.getOutputFileSignatureWithRange());
    evolution.perform(this->out);
}

void Frontend::quench(int argc, char **argv) {
    // Parse options
    cxxopts::Options options(argv[0],
                             Fold("")
                                     .width(80));

    std::string inputFilename;
    std::string outputFilename;
    std::vector<std::string> overridenParamsEntries;
    std::vector<std::string> quenchParamsEntries;
    std::vector<std::string> paramsToPrint{};

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
            ("q,quench_param", "overrides the param as in --set_param, applied after --set_param, but for starting"
                               "Hamiltonian in quantum quench",
             cxxopts::value<std::vector<std::string>>(quenchParamsEntries));

    auto parsedOptions = options.parse(argc, argv);
    if (parsedOptions.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Validate parsed options
    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help");
    if (!parsedOptions.count("input"))
        die("Input file must be specified with option -i [input file name]");
    if (quenchParamsEntries.empty())
        die("At least one parameter should be specified for quantum quench");

    // Prepare quench (initial) parameters
    IO io(std::cout);
    std::vector<std::string> overridenAndQuenchEntries = overridenParamsEntries;
    overridenAndQuenchEntries.insert(overridenAndQuenchEntries.end(), quenchParamsEntries.begin(),
                                     quenchParamsEntries.end());
    Parameters quenchParams = io.loadParameters(inputFilename, overridenAndQuenchEntries);

    // Prepare (final) parameters
    Parameters params = io.loadParameters(inputFilename, overridenParamsEntries);
    for (const auto &param : paramsToPrint)
        if (!params.hasParam(param))
            die("Parameters to print: parameter " + param + " is unknown");

    // Print info on prameters and initial and final Hamiltonians
    params.printGeneral(std::cout);
    std::cout << std::endl;
    std::cout << "[Frontend::quench] Initial Hamiltonian:" << std::endl;
    quenchParams.printHamiltonianTerms(std::cout);
    std::cout << std::endl;
    std::cout << "[Frontend::quench] Final Hamiltonian:" << std::endl;
    params.printHamiltonianTerms(std::cout);
    std::cout << std::endl;

    // Generate Fock basis
    FockBaseGenerator baseGenerator;
    std::cout << "[Frontend::simulate] Preparing Fock basis... " << std::flush;
    arma::wall_clock timer;
    timer.tic();
    auto base = std::shared_ptr(baseGenerator.generate(params.N, params.K));
    std::cout << "done (" << timer.toc() << " s)." << std::endl;

    // Prepare initial and final HamiltonianGenerator
    auto initialRnd = std::make_unique<RND>(params.from + params.seed);
    auto finalRnd = std::make_unique<RND>(params.from + params.seed);
    auto initialHamiltonianGenerator = HamiltonianGeneratorBuilder{}.build(quenchParams, base, *initialRnd);
    auto finalHamiltonianGenerator = HamiltonianGeneratorBuilder{}.build(params, base, *finalRnd);
    auto averagingModel = AveragingModelFactory{}.create(params.averagingModel);

    // Prepare and run quenches
    QuenchCalculator quenchCalculator;
    for (std::size_t i = params.from; i < params.to; i++) {
        std::cout << "[Simulation::quench] Performing quench " << i << "... " << std::flush;
        timer.tic();

        averagingModel->setupHamiltonianGenerator(*initialHamiltonianGenerator, *initialRnd, i, params.totalSimulations);
        averagingModel->setupHamiltonianGenerator(*finalHamiltonianGenerator, *finalRnd, i, params.totalSimulations);
        arma::sp_mat initialHamiltonian = initialHamiltonianGenerator->generate();
        arma::sp_mat finalHamiltonian = finalHamiltonianGenerator->generate();

        quenchCalculator.addQuench(initialHamiltonian, finalHamiltonian);

        std::cout << "done (" << timer.toc() << " s). epsilon: " << quenchCalculator.getLastQuenchEpsilon();
        std::cout << "; std dev: " << quenchCalculator.getLastQuenchEpsilonQuantumUncertainty() << std::endl;
    }

    // Save results
    std::vector<std::string> resultHeader = {"epsilon", "avgError", "quantumError"};
    std::vector<std::string> resultFields = {std::to_string(quenchCalculator.getMeanEpsilon()),
                                             std::to_string(quenchCalculator.getEpsilonAveragingSampleError()),
                                             std::to_string(quenchCalculator.getMeanEpsilonQuantumUncertainty())};

    io.printInlineResults(quenchParams, paramsToPrint, resultHeader, resultFields);
    if (!outputFilename.empty())
        io.storeInlineResults(quenchParams, paramsToPrint, resultHeader, resultFields, outputFilename);
}

void Frontend::printGeneralHelp(const std::string &cmd) {
    std::cout << Fold("Program performing exact diagonalization of Hubbard-like hamiltonians with analyzing "
                      "facilities. ").width(80) << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: " << cmd << " [mode] (mode dependent parameters). " << std::endl;
    std::cout << std::endl;
    std::cout << "Available modes:" << std::endl;
    std::cout << "simulate" << std::endl;
    std::cout << Fold("Performs the diagonalizations according to the parameters passed. Some analyzer tasks can "
                      "be performed on the fly.").width(80).margin(4) << std::endl;
    std::cout << "analyze" << std::endl;
    std::cout << Fold("Performs one or more analyzer tasks after loading simulation results from the files.")
                 .width(80).margin(4) << std::endl;
    std::cout << "chebyshev" << std::endl;
    std::cout << Fold("Performs time evolution using Chebyshev expansion technique.")
                 .width(80).margin(4) << std::endl;
    std::cout << std::endl;
    std::cout << "Type " + cmd + " [mode] --help to get help on the specific mode." << std::endl;
}
