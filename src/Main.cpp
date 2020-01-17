#include <cstdlib>
#include <string>
#include <iostream>
#include <memory>
#include <random>

#include "CavityHamiltonianGenerator.h"
#include "GapRatioCalculator.h"
#include "FockBaseGenerator.h"
#include "Simulation.h"
#include "Utils.h"
#include "Parameters.h"

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

    std::string generate_header(Analyzer &analyzer) {
        std::ostringstream headerStream;
        headerStream << R"("number of sites" "number of bosons" J W U U1 beta phi0 )";

        auto header = analyzer.getInlineResultsHeader();
        std::for_each(header.begin(), header.end(), [](auto &headerField) {
            headerField = "\"" + headerField + "\"";
        });
        std::copy(header.begin(), header.end(), std::ostream_iterator<std::string>(headerStream, " "));
        return headerStream.str();
    }

    std::string generate_fields(const Parameters &params, Analyzer &analyzer) {
        std::ostringstream fieldsStream;
        fieldsStream << params.numberOfSites << " " << params.numberOfBosons << " " << params.J << " ";
        fieldsStream << params.W << " " << params.U << " " << params.U1 << " " << params.beta << " " << params.phi0;
        fieldsStream << " ";

        auto fields = analyzer.getInlineResultsFields();
        std::copy(fields.begin(), fields.end(), std::ostream_iterator<std::string>(fieldsStream, " "));
        return fieldsStream.str();
    }
}

int main(int argc, char **argv) {
    if (argc != 3)
        die(std::string("Usage: ") + argv[0] + " [input file] [output file]");

    std::string inputFilename(argv[1]);
    std::ifstream input(inputFilename);
    if (!input)
        die("Cannot open " + inputFilename + " to read input parameters");

    Parameters params(input);
    params.print(std::cout);
    std::cout << std::endl;

    bool changePhi0ForAverage = (params.phi0 == "changeForAverage");
    auto hamiltonianGenerator = build_hamiltonian_generator(params, changePhi0ForAverage);

    Analyzer analyzer;
    analyzer.addTask(std::make_unique<GapRatioCalculator>(0.5, 0.1));
    if (changePhi0ForAverage) {
        perform_simulations<Phi0AveragingModel>(std::move(hamiltonianGenerator), analyzer, params.numberOfSimulations,
                                                params.saveEigenenergies);
    } else {
        perform_simulations<OnsiteDisorderAveragingModel>(std::move(hamiltonianGenerator), analyzer,
                                                          params.numberOfSimulations, params.saveEigenenergies);
    }

    std::string header = generate_header(analyzer);
    std::string fields = generate_fields(params, analyzer);
    std::cout << std::endl << header << std::endl << fields << std::endl;

    std::string outputFilename(argv[2]);
    save_output_to_file(header, fields, outputFilename);

    return EXIT_SUCCESS;
}


