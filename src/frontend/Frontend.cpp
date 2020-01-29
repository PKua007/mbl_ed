//
// Created by Piotr Kubala on 21/01/2020.
//

#include <random>

#include <cxxopts.hpp>
#include <filesystem>

#include "Frontend.h"
#include "utils/Assertions.h"
#include "simulation/Simulation.h"
#include "simulation/FockBaseGenerator.h"
#include "simulation/CavityHamiltonianGenerator.h"
#include "IO.h"
#include "analyzer/tasks/CDF.h"
#include "analyzer/tasks/MeanInverseParticipationRatio.h"
#include "utils/Fold.h"
#include "utils/Utils.h"

namespace {
    class UniformGenerator {
    private:
        std::mt19937 generator;
        std::uniform_real_distribution<double> distribution;

    public:
        UniformGenerator(double min, double max, unsigned long seed) : generator(seed), distribution(min, max) { }

        UniformGenerator(UniformGenerator &dummy) = delete;
        UniformGenerator operator=(UniformGenerator dummy) = delete;

        double operator()() {
            return this->distribution(this->generator);
        }
    };

    template<typename HamiltonianGenerator_t>
    class OnsiteDisorderAveragingModel {
    public:
        static void setupHamiltonianGenerator(HamiltonianGenerator_t &hamiltonianGenerator, std::size_t simulationIndex,
                                              std::size_t numberOfSimulations)
        {
            static_cast<void>(simulationIndex);
            static_cast<void>(numberOfSimulations);
            hamiltonianGenerator.resampleOnsiteEnergies();
        }
    };

    template<typename HamiltonianGenerator_t>
    class Phi0AveragingModel {
    public:
        static void setupHamiltonianGenerator(HamiltonianGenerator_t &hamiltonianGenerator,
                                              std::size_t simulationIndex, std::size_t numberOfSimulations)
        {
            Expects(numberOfSimulations > 0);
            Expects(simulationIndex < numberOfSimulations);

            double phi0 = 2*M_PI*simulationIndex/numberOfSimulations;
            hamiltonianGenerator.resampleOnsiteEnergies();
            hamiltonianGenerator.setPhi0(phi0);
        }
    };
}

auto Frontend::buildHamiltonianGenerator(const Parameters &params, bool changePhi0ForAverage) {
    FockBaseGenerator baseGenerator;
    auto base = baseGenerator.generate(params.numberOfSites, params.numberOfBosons);

    // We add 'from' to seed not to duplicate results when simulating in parts
    auto disorderGenerator = std::make_unique<UniformGenerator>(-params.W, params.W, params.seed + params.from);

    CavityHamiltonianGeneratorParameters hamiltonianParams;
    hamiltonianParams.J = params.J;
    hamiltonianParams.U = params.U;
    hamiltonianParams.U1 = params.U1;
    hamiltonianParams.beta = params.beta;
    if (!changePhi0ForAverage)
        hamiltonianParams.phi0 = std::stod(params.phi0);

    using TheHamiltonianGenerator = CavityHamiltonianGenerator<UniformGenerator>;
    return std::make_unique<TheHamiltonianGenerator>(std::move(base), hamiltonianParams,
                                                     std::move(disorderGenerator), params.usePeriodicBC);
}

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
        } else if (taskName == "ipr") {
            double mgrCenter, mgrMargin;
            taskStream >> mgrCenter >> mgrMargin;
            ValidateMsg(taskStream, "Wrong format, use: ipr [epsilon center] [epsilon margin]");
            Validate(mgrCenter > 0 && mgrCenter < 1);
            Validate(mgrMargin > 0 && mgrMargin <= 1);
            Validate(mgrCenter - mgrMargin/2 >= 0 && mgrCenter + mgrMargin/2 <= 1);
            analyzer.addTask(std::make_unique<MeanInverseParticipationRatio>(mgrCenter, mgrMargin));
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
void Frontend::perform_simulations(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, Analyzer &analyzer,
                                   const SimulationParameters &simulationParameters)
{
    using TheSimulation = Simulation<HamiltonianGenerator_t, AveragingModel_t<HamiltonianGenerator_t>>;
    TheSimulation simulation(std::move(hamiltonianGenerator), simulationParameters);
    simulation.perform(this->out, analyzer);
}

void Frontend::simulate(int argc, char **argv) {
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
                ->default_value("numberOfSites,numberOfBosons,J,W,U,U1,beta,phi0"));

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help");
    if (!result.count("input"))
        die("Input file must be specified with option -i [input file name]");
    if (!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory))
        die("Output directory " + directory.string() + " does not exist or is not a directory");

    IO io(std::cout);
    Parameters params = io.loadParameters(inputFilename);
    bool changePhi0ForAverage = (params.phi0 == "changeForAverage");
    auto hamiltonianGenerator = buildHamiltonianGenerator(params, changePhi0ForAverage);

    Analyzer analyzer = prepareAnalyzer(onTheFlyTasks);
    std::string fileSignature = params.getOutputFileSignature();
    std::string eigensystemPath = directory / fileSignature;

    SimulationParameters simulationParams;
    simulationParams.from = params.from;
    simulationParams.to = params.to;
    simulationParams.totalSimulations = params.totalSimulations;
    simulationParams.calculateEigenvectors = params.calculateEigenvectors;
    simulationParams.saveEigenenergies = params.saveEigenenergies;
    simulationParams.fileSignature = eigensystemPath;
    if (changePhi0ForAverage)
        perform_simulations<Phi0AveragingModel>(std::move(hamiltonianGenerator), analyzer, simulationParams);
    else
        perform_simulations<OnsiteDisorderAveragingModel>(std::move(hamiltonianGenerator), analyzer, simulationParams);

    io.storeAnalyzerResults(params, analyzer, paramsToPrint, fileSignature, outputFilename);
}

void Frontend::analyze(int argc, char **argv) {
    cxxopts::Options options(argv[0],
                             Fold("Mode for analyzing the results of already performed simulations, stored in the "
                             "files. It performs one or more analyzer task, specified with option -t (--task). The "
                             "files are found using signature generated from parameters.").width(80));

    std::string inputFilename;
    std::string outputFilename;
    std::vector<std::string> paramsToPrint{};
    std::vector<std::string> tasks{};
    std::filesystem::path directory;

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
             cxxopts::value<std::filesystem::path>(directory)->default_value("."));

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    std::string cmd(argv[0]);
    if (argc != 1)
        die("Unexpected positional arguments. See " + cmd + " --help");
    if (!result.count("input"))
        die("Input file must be specified with option -i [input file name]");
    if (!result.count("task"))
        die("At least 1 analyzer task must be specified with option -t [task parameters]");
    if (!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory))
        die("Directory " + directory.string() + " does not exist or is not a directory");

    IO io(std::cout);
    Parameters params = io.loadParameters(inputFilename);
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

    io.storeAnalyzerResults(params, analyzer, paramsToPrint, fileSignature, outputFilename);
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
