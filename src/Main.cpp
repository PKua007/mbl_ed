#include <cstdlib>
#include <string>
#include <iostream>
#include <memory>
#include <random>
#include <filesystem>

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
    Quantity perform_simulations(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                                 std::size_t numberOfSimulations)
    {
        using TheSimulation = Simulation<HamiltonianGenerator_t, AveragingModel_t<HamiltonianGenerator_t>>;
        TheSimulation simulation(std::move(hamiltonianGenerator), numberOfSimulations, 0.5, 0.1);
        simulation.perform(std::cout);
        return simulation.getMeanGapRatio();
    }

    void save_mean_gap_ratio(const Parameters &params, const Quantity &meanGapRatio, const std::string &outputFilename)
    {
        bool fileExists = std::filesystem::exists(outputFilename);
        std::ofstream output(outputFilename, std::fstream::app);
        if (!output)
            die("Cannot open " + outputFilename + " to write mean gap ratio");
        if (!fileExists) {
            output << R"("number of sites" "number of bosons" "J" "W" "U" "U1" "beta" "phi0" "mean gap ratio" )";
            output << R"("mean gap ratio error")" << std::endl;
        }

        output << params.numberOfSites << " " << params.numberOfBosons << " " << params.J << " " << params.W << " ";
        output << params.U << " " << params.U1 << " " << params.beta << " " << params.phi0 << " " << meanGapRatio;
        output << std::endl;
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

    Quantity meanGapRatio;
    if (changePhi0ForAverage) {
        meanGapRatio = perform_simulations<Phi0AveragingModel>(std::move(hamiltonianGenerator),
                                                               params.numberOfSimulations);
    } else {
        meanGapRatio = perform_simulations<OnsiteDisorderAveragingModel>(std::move(hamiltonianGenerator),
                                                                         params.numberOfSimulations);
    }
    meanGapRatio.separator = Quantity::Separator::SPACE;
    std::cout << "Mean gap ratio: " << meanGapRatio << std::endl;

    std::string outputFilename(argv[2]);
    save_mean_gap_ratio(params, meanGapRatio, outputFilename);

    return EXIT_SUCCESS;
}


