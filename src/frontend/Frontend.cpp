//
// Created by Piotr Kubala on 21/01/2020.
//

#include <cxxopts.hpp>
#include <filesystem>
#include <omp.h>

#include "Frontend.h"
#include "IO.h"
#include "CavityConstantsReader.h"

#include "simulation/Simulation.h"
#include "simulation/FockBaseGenerator.h"
#include "simulation/AveragingModels.h"
#include "simulation/DisorderGenerators.h"
#include "simulation/terms/OnsiteDisorder.h"
#include "simulation/terms/HubbardHop.h"
#include "simulation/terms/CavityLongInteraction.h"
#include "simulation/terms/HubbardOnsite.h"
#include "simulation/terms/ListOnsite.h"
#include "simulation/terms/LookupCavityZ2.h"
#include "simulation/terms/LookupCavityYZ.h"
#include "simulation/terms/LookupCavityY2.h"

#include "analyzer/tasks/CDF.h"
#include "analyzer/tasks/MeanInverseParticipationRatio.h"
#include "analyzer/tasks/InverseParticipationRatio.h"
#include "analyzer/tasks/EDCorrelationsTimeEvolution.h"

#include "evolution/ChebyshevEvolution.h"

#include "utils/Fold.h"
#include "utils/Utils.h"
#include "utils/Assertions.h"

/**
 * @brief Builds hamiltonian generator parsing all general parameters and hamiltonian terms from @ params
 */
auto Frontend::buildHamiltonianGenerator(const Parameters &params, std::shared_ptr<FockBase> fockBase, RND &rnd) {
    std::size_t numberOfSites = fockBase->getNumberOfSites();
    auto generator = std::make_unique<HamiltonianGenerator>(fockBase, params.usePeriodicBC);

    for (auto &term : params.hamiltonianTerms) {
        std::string termName = term.first;
        const auto &termParams = term.second;
        if (termName == "hubbardHop") {
            double J = termParams.getDouble("J");
            Validate(J >= 0);
            generator->addHoppingTerm(std::make_unique<HubbardHop>(J));
        } else if (termName == "hubbardOnsite") {
            double U = termParams.getDouble("U");
            Validate(U >= 0);
            generator->addDiagonalTerm(std::make_unique<HubbardOnsite>(U));
        } else if (termName == "onsiteDisorder") {
            double W = termParams.getDouble("W");
            Validate(W >= 0);
            auto disorderGenerator = std::make_unique<UniformGenerator>(-W, W);
            generator->addDiagonalTerm(std::make_unique<OnsiteDisorder<UniformGenerator>>(std::move(disorderGenerator),
                                                                                          numberOfSites, rnd));
        } else if (termName == "listOnsite") {
            auto stringValues = explode(termParams.getString("values"), ',');
            Validate(stringValues.size() == numberOfSites);
            std::vector<double> values(stringValues.size());
            std::transform(stringValues.begin(), stringValues.end(), values.begin(),
                           [](auto s){ return std::stod(s); });
            generator->addDiagonalTerm(std::make_unique<ListOnsite>(values));
        } else if (termName == "quasiperiodicDisorder") {
            double W = termParams.getDouble("W");
            double beta = termParams.getDouble("beta");
            double phi0 = termParams.getDouble("phi0");
            Validate(W >= 0);
            Validate(beta > 0);
            generator->addDiagonalTerm(std::make_unique<QuasiperiodicDisorder>(W, beta, phi0));
        } else if (termName == "cavityLongInteractions") {
            double U1 = termParams.getDouble("U1");
            double beta = termParams.getDouble("beta");
            double phi0 = termParams.getDouble("phi0");
            Validate(U1 >= 0);
            Validate(beta > 0);
            generator->addDiagonalTerm(std::make_unique<CavityLongInteraction>(U1, beta, phi0));
        } else if (termName == "lookupCavityZ2") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addDiagonalTerm(std::make_unique<LookupCavityZ2>(U1, cavityConstants));
        } else if (termName == "lookupCavityYZ") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addHoppingTerm(std::make_unique<LookupCavityYZ>(U1, cavityConstants));
        } else if (termName == "lookupCavityY2") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addDoubleHoppingTerm(std::make_unique<LookupCavityY2>(U1, cavityConstants));
        } else if (termName == "lookupCavityZ2_YZ") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addDiagonalTerm(std::make_unique<LookupCavityZ2>(U1, cavityConstants));
            generator->addHoppingTerm(std::make_unique<LookupCavityYZ>(U1, cavityConstants));
        } else if (termName == "lookupCavityZ2_YZ_Y2") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addDiagonalTerm(std::make_unique<LookupCavityZ2>(U1, cavityConstants));
            generator->addHoppingTerm(std::make_unique<LookupCavityYZ>(U1, cavityConstants));
            generator->addDoubleHoppingTerm(std::make_unique<LookupCavityY2>(U1, cavityConstants));
        } else {
            throw ValidationException("Unknown hamiltonian term: " + termName);
        }
    }

    return generator;
}

/**
 * @brief Takes a vector of @a tasks with parameters and parses them to AnalyzerTask -s
 */
Analyzer Frontend::prepareAnalyzer(const std::vector<std::string> &tasks, const Parameters &params,
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

template<template <typename> typename AveragingModel_t>
void Frontend::perform_simulations(std::unique_ptr<HamiltonianGenerator> hamiltonianGenerator,
                                   std::unique_ptr<RND> rnd, Analyzer &analyzer,
                                   const SimulationParameters &simulationParameters)
{
    using TheSimulation = Simulation<HamiltonianGenerator, AveragingModel_t<UniformGenerator>>;
    TheSimulation simulation(std::move(hamiltonianGenerator), std::move(rnd), simulationParameters);
    simulation.perform(this->out, analyzer);
}

template<template <typename> typename AveragingModel_t>
void Frontend::perform_chebyshev_evolution(std::unique_ptr<HamiltonianGenerator> hamiltonianGenerator,
                                           std::unique_ptr<RND> rnd, const Parameters &params,
                                           const CorrelationsTimeEvolutionParameters &evolutionParameters)
{
    using TheEvolution = ChebyshevEvolution<HamiltonianGenerator, AveragingModel_t<UniformGenerator>>;
    TheEvolution evolution(std::move(hamiltonianGenerator), std::move(rnd), params.from, params.to,
                           params.totalSimulations, evolutionParameters, params.getOutputFileSignatureWithRange());
    evolution.perform(this->out);
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

    FockBaseGenerator baseGenerator;
    auto base = std::shared_ptr(baseGenerator.generate(params.N, params.K));

    auto rnd = std::make_unique<RND>(params.from + params.seed);
    auto hamiltonianGenerator = this->buildHamiltonianGenerator(params, base, *rnd);

    Analyzer analyzer = prepareAnalyzer(onTheFlyTasks, params, base);

    // OpenMP info
    std::cout << "[Frontend::simulate] Using " << omp_get_max_threads() << " OpenMP threads" << std::endl;

    // Prepare and run simulations
    SimulationParameters simulationParams;
    simulationParams.from = params.from;
    simulationParams.to = params.to;
    simulationParams.totalSimulations = params.totalSimulations;
    simulationParams.calculateEigenvectors = params.calculateEigenvectors;
    simulationParams.saveEigenenergies = params.saveEigenenergies;
    simulationParams.fileSignature = directory / params.getOutputFileSignature();
    if (params.averagingModel == "none") {
        perform_simulations<DummyAveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), analyzer,
                                                 simulationParams);
    } else if (params.averagingModel == "uniformPhi0") {
        perform_simulations<UniformPhi0AveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), analyzer,
                                                       simulationParams);
    } else if (params.averagingModel == "randomPhi0") {
        perform_simulations<RandomPhi0AveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), analyzer,
                                                      simulationParams);
    } else if (params.averagingModel == "onsiteDisorder") {
        perform_simulations<OnsiteDisorderAveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), analyzer,
                                                          simulationParams);
    } else if (params.averagingModel == "cavityConstants") {
        perform_simulations<CavityConstantsAveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), analyzer,
                                                           simulationParams);
    } else {
        die("Unknown averaging model: " + params.averagingModel);
    }

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

    FockBaseGenerator baseGenerator;
    auto base = std::shared_ptr(baseGenerator.generate(params.N, params.K));

    // Load eigenenergies and analyze them
    Analyzer analyzer = prepareAnalyzer(tasks, params, base);
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
    io.printInlineAnalyzerResults(params, analyzer, paramsToPrint);
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
    std::string vectorsToEvolveTag;

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
                          " both - both ;)", cxxopts::value<std::string>(vectorsToEvolveTag));

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
        die("You have to specify vectors to evolve using -v [unif/dw/both]");
    // Validation of vectors is done later

    FockBaseGenerator baseGenerator;
    auto base = std::shared_ptr(baseGenerator.generate(params.N, params.K));

    auto rnd = std::make_unique<RND>(params.from + params.seed);
    auto hamiltonianGenerator = this->buildHamiltonianGenerator(params, base, *rnd);

    // OpenMP info
    std::cout << "[Frontend::simulate] Using " << omp_get_max_threads() << " OpenMP threads" << std::endl;

    // Prepare and run evolutions
    CorrelationsTimeEvolutionParameters evolutionParameters;

    std::istringstream timeSegmentationStream(timeSegmentation);
    std::copy(std::istream_iterator<EvolutionTimeSegment>(timeSegmentationStream),
              std::istream_iterator<EvolutionTimeSegment>(),
              std::back_inserter(evolutionParameters.timeSegmentation));
    evolutionParameters.numberOfSites = params.K;
    evolutionParameters.fockBase = base;
    evolutionParameters.marginSize = marginSize;
    evolutionParameters.setVectorsToEvolveFromTag(vectorsToEvolveTag); // This one also does the validation
    if (params.averagingModel == "none") {
        perform_chebyshev_evolution<DummyAveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), params,
                                                         evolutionParameters);
    } else if (params.averagingModel == "uniformPhi0") {
        perform_chebyshev_evolution<UniformPhi0AveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), params,
                                                               evolutionParameters);
    } else if (params.averagingModel == "randomPhi0") {
        perform_chebyshev_evolution<RandomPhi0AveragingModel>(std::move(hamiltonianGenerator), std::move(rnd), params,
                                                              evolutionParameters);
    } else if (params.averagingModel == "onsiteDisorder") {
        perform_chebyshev_evolution<OnsiteDisorderAveragingModel>(std::move(hamiltonianGenerator), std::move(rnd),
                                                                  params, evolutionParameters);
    } else if (params.averagingModel == "cavityConstants") {
        perform_chebyshev_evolution<CavityConstantsAveragingModel>(std::move(hamiltonianGenerator), std::move(rnd),
                                                                   params, evolutionParameters);
    } else {
        die("Unknown averaging model: " + params.averagingModel);
    }
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
