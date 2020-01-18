#include <cstdlib>
#include <string>
#include <iostream>
#include <memory>
#include <random>

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

    Analyzer prepare_analyzer(const std::string &onTheFlyTasks) {
        auto tasks = explode(onTheFlyTasks, ';');
        std::for_each(tasks.begin(), tasks.end(), trim);
        tasks.erase(std::remove_if(tasks.begin(), tasks.end(), std::mem_fn(&std::string::empty)), tasks.end());

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
}

int main(int argc, char **argv) {
    if (argc != 4 && argc != 5) {
        die(std::string("Usage: ") + argv[0] + " [input file] [on the fly tasks] [output file] "
            + "(parameters to print = all)");
    }

    std::string inputFilename(argv[1]);
    std::ifstream input(inputFilename);
    if (!input)
        die("Cannot open " + inputFilename + " to read input parameters");

    Parameters params(input);
    params.print(std::cout);
    std::cout << std::endl;

    bool changePhi0ForAverage = (params.phi0 == "changeForAverage");
    auto hamiltonianGenerator = build_hamiltonian_generator(params, changePhi0ForAverage);

    std::string onTheFlyTasks(argv[2]);

    Analyzer analyzer = prepare_analyzer(onTheFlyTasks);
    std::string fileSignature = hamiltonianGenerator->fileSignature();
    if (changePhi0ForAverage) {
        perform_simulations<Phi0AveragingModel>(std::move(hamiltonianGenerator), analyzer, params.numberOfSimulations,
                                                params.saveEigenenergies);
    } else {
        perform_simulations<OnsiteDisorderAveragingModel>(std::move(hamiltonianGenerator), analyzer,
                                                          params.numberOfSimulations, params.saveEigenenergies);
    }

    std::vector<std::string> paramsToPrint;
    if (argc == 5) {
        paramsToPrint = explode(argv[4], ' ');
        paramsToPrint.erase(std::remove_if(paramsToPrint.begin(), paramsToPrint.end(),
                                           std::mem_fn(&std::string::empty)),
                            paramsToPrint.end());
    } else {
        paramsToPrint = {"numberOfSites", "numberOfBosons", "J", "W", "U", "U1", "beta", "phi0"};
    }
    InlineResultsPrinter resultsPrinter(params, analyzer, paramsToPrint);
    std::cout << std::endl << resultsPrinter.getHeader() << std::endl << resultsPrinter.getFields() << std::endl;
    std::string outputFilename(argv[3]);
    save_output_to_file(resultsPrinter.getHeader(), resultsPrinter.getFields(), outputFilename);

    std::cout << std::endl << "Storing bulk results... " << std::flush;
    analyzer.storeBulkResults(fileSignature);
    std::cout << "done." << std::endl;

    return EXIT_SUCCESS;
}


