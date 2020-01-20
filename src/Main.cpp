#include <cstdlib>
#include <string>
#include <iostream>
#include <memory>
#include <random>
#include <filesystem>

#include <cxxopts.hpp>

#include "simulation/CavityHamiltonianGenerator.h"
#include "analyzer/tasks/MeanGapRatio.h"
#include "analyzer/tasks/CDF.h"
#include "simulation/FockBaseGenerator.h"
#include "simulation/Simulation.h"
#include "utils/Utils.h"
#include "Parameters.h"
#include "InlineResultsPrinter.h"

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

    auto build_hamiltonian_generator(const Parameters &params, bool changePhi0ForAverage) {
        FockBaseGenerator baseGenerator;
        auto base = baseGenerator.generate(params.numberOfSites, params.numberOfBosons);

        auto disorderGenerator = std::make_unique<UniformGenerator>(-params.W, params.W, params.seed);

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

    template<template<typename> typename AveragingModel_t, typename HamiltonianGenerator_t>
    void perform_simulations(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, Analyzer &analyzer,
                             std::size_t numberOfSimulations, bool saveEigenenergies)
    {
        using TheSimulation = Simulation<HamiltonianGenerator_t, AveragingModel_t<HamiltonianGenerator_t>>;
        TheSimulation simulation(std::move(hamiltonianGenerator), numberOfSimulations, saveEigenenergies);
        simulation.perform(std::cout, analyzer);
    }

    void save_output_to_file(const std::string &header, const std::string &fields, const std::string &outputFilename) {
        bool fileExists = std::ifstream(outputFilename).is_open();
        std::ofstream output(outputFilename, std::fstream::app);
        if (!output)
            die("Cannot open " + outputFilename + " to write analyzer output");
        if (!fileExists)
            output << header << std::endl;

        output << fields << std::endl;
    }

    Analyzer prepare_analyzer(const std::vector<std::string> &tasks) {
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

    void store_analyzer_results(const Parameters &params, const Analyzer &analyzer,
                                const std::vector<std::string> &paramsToPrint, const std::string &fileSignature,
                                const std::string &outputFilename)
    {
        InlineResultsPrinter resultsPrinter(params, analyzer.getInlineResultsHeader(),
                                            analyzer.getInlineResultsFields(), paramsToPrint);
        std::cout << std::endl << resultsPrinter.getHeader() << std::endl << resultsPrinter.getFields() << std::endl;

        if (!outputFilename.empty())
            save_output_to_file(resultsPrinter.getHeader(), resultsPrinter.getFields(), outputFilename);

        std::cout << std::endl << "Storing bulk results... " << std::flush;
        analyzer.storeBulkResults(fileSignature);
        std::cout << "done." << std::endl;
    }

    Parameters load_parameters(const std::string &inputFilename) {
        std::ifstream input(inputFilename);
        if (!input)
            die("Cannot open " + inputFilename + " to read input parameters");

        Parameters params(input);
        params.print(std::cout);
        std::cout << std::endl;
        return params;
    }

    void simulate(int argc, char **argv) {
        cxxopts::Options options(argv[0], " - exact diagonalization mode");

        std::string inputFilename;
        std::string outputFilename;
        std::vector<std::string> paramsToPrint;
        std::vector<std::string> onTheFlyTasks;

        options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters", cxxopts::value<std::string>(inputFilename))
            ("o,output", "when specified, inline results will be printed to this file",
                cxxopts::value<std::string>(outputFilename))
            ("t,task", "task(s) to be performed on the fly while simulating",
                cxxopts::value<std::vector<std::string>>(onTheFlyTasks))
            ("p,print_parameter", "parameters to be included in inline results",
                cxxopts::value<std::vector<std::string>>(paramsToPrint)
                    ->default_value("numberOfSites,numberOfBosons,J,W,U,U1,beta,phi0"));

        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        if (!result.count("input"))
            die("Input file must be specified with option -i [input file name]");

        Parameters params = load_parameters(inputFilename);
        bool changePhi0ForAverage = (params.phi0 == "changeForAverage");
        auto hamiltonianGenerator = build_hamiltonian_generator(params, changePhi0ForAverage);

        Analyzer analyzer = prepare_analyzer(onTheFlyTasks);
        std::string fileSignature = hamiltonianGenerator->fileSignature();
        if (changePhi0ForAverage) {
            perform_simulations<Phi0AveragingModel>(std::move(hamiltonianGenerator), analyzer,
                                                    params.numberOfSimulations, params.saveEigenenergies);
        } else {
            perform_simulations<OnsiteDisorderAveragingModel>(std::move(hamiltonianGenerator), analyzer,
                                                              params.numberOfSimulations, params.saveEigenenergies);
        }

        store_analyzer_results(params, analyzer, paramsToPrint, fileSignature, outputFilename);
    }

    std::vector<std::string> find_eigenenergy_files(const std::string &directory, const std::string &fileSignature) {
        std::vector<std::string> files;
        for (const auto &entry : std::filesystem::directory_iterator(directory)) {
            std::string filename = entry.path();
            std::string prefix = directory + "/" + fileSignature + "_";
            if (startsWith(filename, prefix) && endsWith(filename, "_nrg.dat"))
                files.push_back(filename);
        }

        return files;
    }

    void analyze(int argc, char **argv) {
        cxxopts::Options options(argv[0], " - analyze mode");

        std::string inputFilename;
        std::string outputFilename;
        std::vector<std::string> paramsToPrint;
        std::vector<std::string> tasks;
        std::string directory;
        std::string fileSignature;

        options.add_options()
            ("h,help", "prints help for this mode")
            ("i,input", "file with parameters", cxxopts::value<std::string>(inputFilename))
            ("o,output", "when specified, inline results will be printed to this file",
                cxxopts::value<std::string>(outputFilename))
            ("t,task", "analyzer task(s) to be performed",
                cxxopts::value<std::vector<std::string>>(tasks))
            ("p,print_parameter", "parameters to be included in inline results",
                cxxopts::value<std::vector<std::string>>(paramsToPrint)
                     ->default_value("numberOfSites,numberOfBosons,J,W,U,U1,beta,phi0"))
            ("d,directory", "directory to search simulation results",
                cxxopts::value<std::string>(directory)->default_value("."))
            ("f,file_signature", "the signature of result files",
                 cxxopts::value<std::string>(fileSignature));

        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }
        
        if (!result.count("input"))
            die("Input file must be specified with option -i [input file name]");
        if (!result.count("task"))
            die("At least 1 analyzer task must be specified with option -t [task parameters]");
        if (!result.count("file_signature"))
            die("File signature must be specified with option -f [file signature]");

        Parameters params = load_parameters(inputFilename);
        Analyzer analyzer = prepare_analyzer(tasks);
        std::vector<std::string> energiesFilenames = find_eigenenergy_files(directory, fileSignature);
        if (energiesFilenames.empty())
            die("No eigenenergy files were found.");

        for (const auto &energiesFilename : energiesFilenames) {
            std::ifstream energiesFile(energiesFilename);
            if (!energiesFile)
                die("Cannot open " + energiesFilename + " to read eigenenergies from");

            std::vector<double> eigenenergies;
            std::copy(std::istream_iterator<double>(energiesFile), std::istream_iterator<double>(),
                      std::back_inserter(eigenenergies));

            analyzer.analyze(eigenenergies);
        }

        store_analyzer_results(params, analyzer, paramsToPrint, fileSignature, outputFilename);
    }
}

int main(int argc, char **argv) {
    std::string cmd(argv[0]);
    if (argc < 2)
        die("Usage: " + cmd + " [mode] (mode dependent parameters)");

    // We now shift the arguments, and pretend, that the new first (command) is "cmd mode"
    // Different modes can then parse the arguments separately
    std::string mode(argv[1]);
    std::string cmdAndMode = cmd + " " + mode;
    argv++;
    argc--;
    argv[0] = cmdAndMode.data();

    if (mode == "simulate") {
        simulate(argc, argv);
    } else if (mode == "analyze") {
        analyze(argc, argv);
    } else {
        die("Unknown mode " + mode);
    }

    return EXIT_SUCCESS;
}




