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

        using TheHamiltonianGenerator = CavityHamiltonianGenerator<UniformGenerator>;
        auto disorderGenerator = std::make_unique<UniformGenerator>(-params.W, params.W, params.seed);

        CavityHamiltonianGeneratorParameters hamiltonianParams;
        hamiltonianParams.J = params.J;
        hamiltonianParams.U = params.U;
        hamiltonianParams.U1 = params.U1;
        hamiltonianParams.beta = params.beta;
        if (!changePhi0ForAverage)
            hamiltonianParams.phi0 = std::stod(params.phi0);

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
}

int main(int argc, char **argv) {
    if (argc != 2)
        die(std::string("Usage: ") + argv[0] + " [input file]");

    std::string filename(argv[1]);
    std::ifstream input(filename);
    if (!input)
        die("Cannot open " + filename + " to read input parameters");

    Parameters params(input);
    params.print(std::cout);
    std::cout << std::endl;

    bool changePhi0ForAverage = (params.phi0 == "changeForAverage");
    auto hamiltonianGenerator = build_hamiltonian_generator(params, changePhi0ForAverage);

    Quantity meanGapRatio;
    meanGapRatio.separator = Quantity::Separator::PLUS_MINUS;
    if (changePhi0ForAverage) {
        meanGapRatio = perform_simulations<Phi0AveragingModel>(std::move(hamiltonianGenerator),
                                                               params.numberOfSimulations);
    } else {
        meanGapRatio = perform_simulations<OnsiteDisorderAveragingModel>(std::move(hamiltonianGenerator),
                                                                         params.numberOfSimulations);
    }
    std::cout << "Mean gap ratio: " << meanGapRatio << std::endl;

    return EXIT_SUCCESS;
}