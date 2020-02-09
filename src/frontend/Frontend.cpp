//
// Created by Piotr Kubala on 21/01/2020.
//

#include <cxxopts.hpp>
#include <filesystem>

#include "Frontend.h"
#include "IO.h"

#include "simulation/Simulation.h"
#include "simulation/FockBaseGenerator.h"
//#include "simulation/CavityHamiltonianGenerator.h"
#include "simulation/AveragingModels.h"
#include "simulation/DisorderGenerators.h"

#include "analyzer/tasks/CDF.h"
#include "analyzer/tasks/MeanInverseParticipationRatio.h"
#include "analyzer/tasks/InverseParticipationRatio.h"

#include "utils/Fold.h"
#include "utils/Utils.h"
#include "utils/Assertions.h"

/**
 * @param changePhi0ForAverage this determines whether to average on disorder or phi0
 */
/*auto Frontend::buildHamiltonianGenerator(const Parameters &params, RND &rnd) {
    FockBaseGenerator baseGenerator;
    auto base = baseGenerator.generate(params.numberOfSites, params.numberOfBosons);

    // We add 'from' to seed not to duplicate results when simulating in parts
    auto disorderGenerator = std::make_unique<UniformGenerator>(-params.W, params.W);

    CavityHamiltonianParameters hamiltonianParams;
    hamiltonianParams.J = params.J;
    hamiltonianParams.U = params.U;
    hamiltonianParams.U1 = params.U1;
    hamiltonianParams.beta = params.beta;
    hamiltonianParams.phi0 = params.phi0;

    using TheHamiltonianGenerator = CavityHamiltonianGenerator<UniformGenerator>;
    return std::make_unique<TheHamiltonianGenerator>(std::move(base), hamiltonianParams, rnd,
                                                     std::move(disorderGenerator), params.usePeriodicBC);
}*/

/**
 * @brief Takes a vector of @a tasks with parameters and parses them to AnalyzerTask -s
 */
Analyzer Frontend::prepareAnalyzer(const std::vector<std::string> &tasks) {
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
        } else {
            throw ValidationException("Unknown analyzer task: " + taskName);
        }
    }
    return analyzer;
}

template<template <typename> typename AveragingModel_t, typename HamiltonianGenerator_t>
void Frontend::perform_simulations(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                                   std::unique_ptr<RND> rnd, Analyzer &analyzer,
                                   const SimulationParameters &simulationParameters)
{
    using TheSimulation = Simulation<HamiltonianGenerator_t, AveragingModel_t<HamiltonianGenerator_t>>;
    TheSimulation simulation(std::move(hamiltonianGenerator), std::move(rnd), simulationParameters);
    simulation.perform(this->out, analyzer);
}

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
            ("i,input", "file with parameters. See input.txt for parameters description",
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
                ->default_value("numberOfSites,numberOfBosons,J,W,U,U1,beta,phi0"))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P J=1 (-PJ=1 does not work) act as one would append J=1 to input file",
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

    auto rnd = std::make_unique<RND>(params.from + params.seed);
    //auto hamiltonianGenerator = this->buildHamiltonianGenerator(params, *rnd);

    Analyzer analyzer = prepareAnalyzer(onTheFlyTasks);

    // Prepare and run simulations
    SimulationParameters simulationParams;
    simulationParams.from = params.from;
    simulationParams.to = params.to;
    simulationParams.totalSimulations = params.totalSimulations;
    simulationParams.calculateEigenvectors = params.calculateEigenvectors;
    simulationParams.saveEigenenergies = params.saveEigenenergies;
    simulationParams.fileSignature = directory / params.getOutputFileSignature();
    /*if (params.averagingModel == "uniformPhi0") {
        perform_simulations<UniformPhi0AveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), analyzer,
                                                       simulationParams);
    } else if (params.averagingModel == "randomPhi0") {
        perform_simulations<RandomPhi0AveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), analyzer,
                                                      simulationParams);
    } else if (params.averagingModel == "onsiteDisorder") {
        perform_simulations<OnsiteDisorderAveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), analyzer,
                                                          simulationParams);
    }*/

    // Save results
    io.printInlineAnalyzerResults(params, analyzer, paramsToPrint);
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
            ("i,input", "file with parameters. See input.txt for parameters description",
             cxxopts::value<std::string>(inputFilename))
            ("o,output", "when specified, inline results will be printed to this file",
             cxxopts::value<std::string>(outputFilename))
            ("t,task", "analyzer task(s) to be performed",
             cxxopts::value<std::vector<std::string>>(tasks))
            ("p,print_parameter", "parameters to be included in inline results",
             cxxopts::value<std::vector<std::string>>(paramsToPrint)
                     ->default_value("numberOfSites,numberOfBosons,J,W,U,U1,beta,phi0"))
            ("d,directory", "directory to search simulation results",
             cxxopts::value<std::filesystem::path>(directory)->default_value("."))
            ("P,set_param", "overrides the value of the parameter loaded as --input. More precisely, doing "
                            "-P J=1 (-PJ=1 does not work) act as one would append J=1 to input file",
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

    // Load eigenenergies and analyze them
    Analyzer analyzer = prepareAnalyzer(tasks);
    std::string fileSignature = params.getOutputFileSignature();
    std::vector<std::string> energiesFilenames = io.findEigenenergyFiles(directory, fileSignature);
    if (energiesFilenames.empty())
        die("No eigenenergy files were found.");

    for (const auto &energiesFilename : energiesFilenames) {
        std::ifstream energiesFile(energiesFilename);
        if (!energiesFile)
            die("Cannot open " + energiesFilename + " to read eigenenergies from");

        Eigensystem eigensystem;
        eigensystem.restore(energiesFile);
        analyzer.analyze(eigensystem);
    }

    // Save results
    io.printInlineAnalyzerResults(params, analyzer, paramsToPrint);
    if (outputFilename.empty())
        io.storeAnalyzerResults(params, analyzer, paramsToPrint, std::nullopt);
    else
        io.storeAnalyzerResults(params, analyzer, paramsToPrint, outputFilename);
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
    std::cout << std::endl;
    std::cout << "Type " + cmd + " [mode] --help to get help on the specific mode." << std::endl;
}
